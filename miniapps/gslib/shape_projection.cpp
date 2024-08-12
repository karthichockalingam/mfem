// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.
//
//      --------------------------------------------------------------
//      Field Diff Miniapp: Compare grid functions on different meshes
//      --------------------------------------------------------------
//
// This miniapp compares two different high-order grid functions, defined on two
// different high-order meshes, based on the GSLIB-FindPoints general off-grid
// interpolation utility. Using a set of points defined within the bounding box
// of the domain, FindPoints is used to interpolate the grid functions from the
// two different meshes and output the difference between the interpolated
// values. The miniapp also uses FindPoints to interpolate the solution from one
// mesh onto another, and visualize the difference using GLVis.
//
// Compile with: make field-diff
//
// Sample runs:
//    field-diff
//    field-diff -m1 triple-pt-1.mesh -s1 triple-pt-1.gf -m2 triple-pt-2.mesh -s2 triple-pt-2.gf -p 200

#include "mfem.hpp"
#include <fstream>

using namespace mfem;
using namespace std;

int main (int argc, char *argv[])
{

   Mesh mesh_macro = Mesh::MakeCartesian2D(2, 2, mfem::Element::Type::QUADRILATERAL,
                                     true, 2.0, 2.0);
   
   Mesh mesh_micro = Mesh::MakeCartesian2D(4, 4, mfem::Element::Type::QUADRILATERAL,
                                     true, 0.8, 0.8, true, 0.1, 0.1);

   int dim = mesh_macro.Dimension();
/*
   Vector vert_coord;

   mesh_micro.GetNodes(vert_coord);
   int nv = vert_coord.Size()/dim;

   for (int i = 0; i < nv; i++)
   {
      vert_coord(0*nv+i) += 0.1;
      vert_coord(1*nv+i) += 0.1;
   }
*/
   Vector vert_coord_copy;

   mesh_micro.GetNodes(vert_coord_copy);
   int nv = vert_coord_copy.Size()/dim;

   for (int i = 0; i < nv; i++)
   {
      std::cout << i << "  " << vert_coord_copy(0*nv+i) << "  " << vert_coord_copy(1*nv+i) << std::endl;
   }

   int order = 1;
   mesh_micro.SetCurvature(order, false, dim, 0); 
   mesh_macro.SetCurvature(order, false, dim, 0); 

   FiniteElementCollection *fec = new H1_FECollection(1, dim);

   FiniteElementSpace fespace_macro(&mesh_macro, fec);
   FiniteElementSpace fespace_micro(&mesh_micro, fec);

   Array<int> boundary_dofs;
   fespace_micro.GetBoundaryTrueDofs(boundary_dofs);

   Array<int> dofs; //to store edge dofs

   set<int> dof_set;
   
   int bnd_pts =  boundary_dofs.Size();
   Vector vxy(bnd_pts * dim);
   int* dof_ids = new int[bnd_pts];
   int count = 0;

   for (int i = 0; i < mesh_micro.GetNBE(); i++)
   {
      const FiniteElement * elem = fespace_micro.GetBE(i);
      const IntegrationRule & ir = elem->GetNodes();
       
      fespace_micro.GetBdrElementDofs(i, dofs);

      FaceElementTransformations * T = mesh_micro.GetBdrFaceTransformations(i);
      DenseMatrix P;
      T->Transform(ir, P);

      assert(("The size of dof should be same ir", dofs.Size()  == ir.Size()));
      
      for (int j = 0; j < ir.Size(); j++)
      {
       // std::cout << "DOF# " << dofs[j] << "  "  << P(0, j) << "  " << P(1, j) << std::endl;
  
       if(dof_set.find(dofs[j]) == dof_set.end())
        {
        vxy(count)             = P(0, j);
        vxy(bnd_pts + (count)) = P(1, j);
        dof_set.insert(dofs[j]); 
        dof_ids[count] = dofs[j];
        count = count + 1;
        }
      }
   }
  
   FindPointsGSLIB finder_pos;
   finder_pos.FindPoints(mesh_macro, vxy);
   
   auto obj_elem = finder_pos.GetGSLIBElem();
   auto obj_ref_pos = finder_pos.GetReferencePosition();
   //auto obj_ref_pos = finder_pos.GetGSLIBReferencePosition();

   //std::cout << "Boundary points:" << bnd_pts << " Count = " << count << std::endl;

   IntegrationPoint ip;
   BiLinear2DFiniteElement bilinear_elem;
   Vector shape(dim*dim);

   GridFunction u_micro_gf(&fespace_micro);
   u_micro_gf = 0.0;

   for (int i = 0; i < bnd_pts; i++)
   {
      ip.Set2(obj_ref_pos(i*dim + 0), obj_ref_pos(i*dim + 1));
      bilinear_elem.CalcShape(ip, shape);

      //std::cout << obj_elem[i] << " DOF # " << dof_ids[i]  << "  " << obj_ref_pos(i*dim + 0) << "  " << obj_ref_pos(i*dim + 1) << "  " << shape(0) << std::endl;
      std::cout << " DOF # " << dof_ids[i]  << "  x =  " <<  vxy(i) << "  y =   " << vxy(bnd_pts + i) << "  macro_value =  " << shape(0) << std::endl;
      u_micro_gf(dof_ids[i]) = shape(0);
   }

   GridFunctionCoefficient u_micro_cgf(&u_micro_gf);

   Array<int> ess_tdof_list;
   Array<int> ess_bdr(mesh_micro.bdr_attributes.Max());
   ess_bdr = 1;
   fespace_micro.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);


   GridFunction x(&fespace_micro);
   x = 0.0;

   x.ProjectBdrCoefficient(u_micro_cgf, ess_bdr);

   // 6. Set up the linear form b(.) corresponding to the right-hand side.
   ConstantCoefficient zero(0.0);
   LinearForm b(&fespace_micro);
   b.AddDomainIntegrator(new DomainLFIntegrator(zero));
   b.Assemble();

   // 7. Set up the bilinear form a(.,.) corresponding to the -Delta operator.
   BilinearForm a(&fespace_micro);
   a.AddDomainIntegrator(new DiffusionIntegrator);
   a.Assemble();

   // 8. Form the linear system A X = B. This includes eliminating boundary
   //    conditions, applying AMR constraints, and other transformations.
   SparseMatrix A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

   // 9. Solve the system using PCG with symmetric Gauss-Seidel preconditioner.
   GSSmoother M(A);
   PCG(A, M, B, X, 1, 200, 1e-12, 0.0);

   // 10. Recover the solution x as a grid function and save to file. The output
   //     can be viewed using GLVis as follows: "glvis -m mesh.mesh -g sol.gf"
   a.RecoverFEMSolution(X, b, x);

   ParaViewDataCollection pd("Diffusion", &mesh_micro);
   pd.SetPrefixPath("ParaView");
   pd.RegisterField("solution", &x);
   pd.SetLevelsOfDetail(order);
   pd.SetDataFormat(VTKFormat::BINARY);
   pd.SetHighOrderOutput(true);
   pd.SetCycle(0);
   pd.SetTime(0.0);
   pd.Save();
   
   return 0;
}
