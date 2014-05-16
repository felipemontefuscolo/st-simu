#include "common.hpp"

PetscErrorCode FormJacobian(SNES snes, Vec Vec_u, Mat *Mat_jac, Mat *prejac, MatStructure *flag, void *ptr);
PetscErrorCode FormFunction(SNES snes, Vec Vec_u, Vec Vec_res, void *ptr);
PetscErrorCode CheckSnesConvergence(SNES snes, PetscInt it, PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx);


void AppCtx::initAll()
{
  initMesh();
  initShapes();
  initMeshIo();
  initDofMappers();
  initPetscObjs();
}

void AppCtx::initMesh()
{
  if (mp)
    delete mp;

  mp = new MeshT(settings.space_dim);
  
  mesh_reader.readFile(settings.input_meshfile.c_str(), mp);
  //

}

void AppCtx::initMeshIo()
{
  mesh_printer.attachMesh(mp);
}

void AppCtx::initShapes()
{
  shapes_c.resize(unk_fields.size());
  shapes_f.resize(unk_fields.size());
  shapes_r.resize(unk_fields.size());

  // defining the shape functions
  for (unsigned i = 0; i < shapes_c.size(); ++i)
  {
    shapes_c[i].setType( unk_fields[i].interpolation.c_str(), cell_dim(CELL_TYPE),   unk_fields[i].degree);
    shapes_f[i].setType( unk_fields[i].interpolation.c_str(), cell_dim(FACET_TYPE),  unk_fields[i].degree);
    if (cell_dim(CELL_TYPE) > 2)
      shapes_r[i].setType( unk_fields[i].interpolation.c_str(), cell_dim(RIDGE_TYPE),  unk_fields[i].degree);
  }

  // defining the quadrature points
  Q_c.setType(CELL_TYPE, cell_dim(CELL_TYPE));
  Q_f.setType(FACET_TYPE, cell_dim(FACET_TYPE));
  if (cell_dim(CELL_TYPE) > 2)
    Q_r.setType(RIDGE_TYPE, cell_dim(RIDGE_TYPE));

  // allocating
  phi_c.reshape(Q_c.numPoints(), shapes_c.size());
  phi_f.reshape(Q_f.numPoints(), shapes_f.size());
  if (cell_dim(CELL_TYPE) > 2)
    phi_r.reshape(Q_r.numPoints(), shapes_r.size());

  // allocating
  dLphi_c.reshape(Q_c.numPoints(), shapes_c.size());
  dLphi_f.reshape(Q_f.numPoints(), shapes_f.size());
  if (cell_dim(CELL_TYPE) > 2)
    dLphi_r.reshape(Q_r.numPoints(), shapes_r.size());


  // evaluating shape-functions at quadrature points
  // cell
  for (unsigned qp = 0; qp < phi_c.dim(0); ++qp)
    for (unsigned var = 0; var < phi_c.dim(1); ++var)
    {
      phi_c[qp][var].resize(shapes_c[var].numDofs());
      dLphi_c[qp][var].reshape(shapes_c[var].numDofs(), cell_dim(CELL_TYPE));

      for (int i = 0; i < shapes_c[var].numDofs(); ++i)
      {
        double const* x = Q_c.point(qp);
        phi_c[qp][var][i] = shapes_c[var].value(x, i);
        for (int c = 0; c < cell_dim(CELL_TYPE); ++c)
          dLphi_c[qp][var][i][c] =  shapes_c[var].grad(x, i, c);
      }
    }

  // facet
  for (unsigned qp = 0; qp < phi_f.dim(0); ++qp)
    for (unsigned var = 0; var < phi_f.dim(1); ++var)
    {
      phi_f[qp][var].resize(shapes_f[var].numDofs());
      dLphi_f[qp][var].reshape(shapes_f[var].numDofs(), cell_dim(FACET_TYPE));

      for (int i = 0; i < shapes_f[var].numDofs(); ++i)
      {
        double const* x = Q_f.point(qp);
        phi_f[qp][var][i] = shapes_f[var].value(x, i);
        for (int c = 0; c < cell_dim(FACET_TYPE); ++c)
          dLphi_f[qp][var][i][c] =  shapes_f[var].grad(x, i, c);
      }
    }

  // ridge
  if (cell_dim(FACET_TYPE) > 2)
    for (unsigned qp = 0; qp < phi_r.dim(0); ++qp)
      for (unsigned var = 0; var < phi_r.dim(1); ++var)
      {
        phi_r[qp][var].resize(shapes_r[var].numDofs());
        dLphi_r[qp][var].reshape(shapes_r[var].numDofs(), cell_dim(RIDGE_TYPE));

        for (int i = 0; i < shapes_r[var].numDofs(); ++i)
        {
          double const* x = Q_r.point(qp);
          phi_r[qp][var][i] = shapes_r[var].value(x, i);
          for (int c = 0; c < cell_dim(RIDGE_TYPE); ++c)
            dLphi_r[qp][var][i][c] =  shapes_r[var].grad(x, i, c);
        }
      }



}

void AppCtx::initDofMappers()
{
  if (!mp) {
    printf("ERROR: Must create the mesh first!\n");
    throw;
  }

  // assigning these system to the correspondent variables
  for (unsigned i = 0; i < systems.size(); ++i)
    for (unsigned j = 0; j < systems[i].fields_id.size(); ++j)
      unk_fields[systems[i].fields_id[j]].sys_num = i;

  
  for (unsigned i = 0; i < systems.size(); ++i)
  {
    if (systems[i].fields_id.empty())
    {
      printf("ERROR: can not create a system with no fields (variables)\n");
      throw;
    }
    
    // create its own mapper
    if (systems[i].sys_mapper < 0)
    {
      systems[i].mapper.reset(new DofMapperT(mp));
      
      for (unsigned j = 0; j < systems[i].fields_id.size(); ++j)
      {
        int const vid = systems[i].fields_id[j];
        AppCtx::UnkField& v = unk_fields[vid];
        systems[i].mapper->addVariable(v.name.c_str(), v.n_comps*shapes_c[vid].numDofsPerVertex(),
                                                       v.n_comps*shapes_c[vid].numDofsInRidge(),
                                                       v.n_comps*shapes_c[vid].numDofsInFacet(),
                                                       v.n_comps*shapes_c[vid].numDofsInCell(),
                                                       v.n_regions, v.n_tags.data(), v.tags.data());
      }
      
      systems[i].mapper->SetUp(systems[i].first_dof_number);
      systems[i].n_dofs = systems[i].mapper->numDofs() + systems[i].first_dof_number;
    }
  }
  
  // Reuse mappers for systems specified by the user
  for (unsigned i = 0; i < systems.size(); ++i)
  {
    int other_sys = systems[i].sys_mapper;
    
    if (other_sys < 0)
      continue;
    
    if (other_sys >= (int)systems.size() || (other_sys==(int)i))
    {
      printf("ERROR: in input file: in system %d: invalid system ID (mapper)\n", (int)i);
      throw;
    }
    
    systems[i].mapper = systems[other_sys].mapper;
    systems[i].n_dofs = systems[i].mapper->numDofs() + systems[i].first_dof_number;
  }

}


// It is only a estimative ... it could be better if the nnz is exact
int estimateNnz(MeshT const* mp, DofMapperT const& mapper)
{
  // get the greater valency
  std::vector<CellH> star, star2;
  VertexH v = mp->vertexBegin();
  VertexH const v_end = mp->vertexEnd();
  for (; v != v_end; ++v)
  {
    if (v.isDisabled(mp))
      continue;

    star2 = v.star(mp);
    if (star.size() < star2.size())
      star = star2;
  }

  int nnz = 0;
  // vertices
  {
    SetVector<int> x;
    for (unsigned i = 0; i < star.size(); ++i)
    {
      for (int j = 0; j < MeshT::verts_per_cell; ++j)
        x.insert( star[i].vertex(mp, j).id(mp) );
    }
    // # dofs
    for (int i = 0; i < mapper.numVars(); ++i)
    {
      nnz += x.size()*mapper.variable(i).numDofsInVertex();
    }
  }
  
  // ridges
  if (MeshT::cell_dim < 2)
  {
    SetVector<int> x;
    for (unsigned i = 0; i < star.size(); ++i)
    {
      for (int j = 0; j < MeshT::ridges_per_cell; ++j)
        x.insert( star[i].ridge(mp, j).id(mp) );
    }
    // # dofs
    for (int i = 0; i < mapper.numVars(); ++i)
    {
      nnz += x.size()*mapper.variable(i).numDofsInRidge();
    }
  }
  
  // facets
  {
    SetVector<int> x;
    for (unsigned i = 0; i < star.size(); ++i)
    {
      for (int j = 0; j < MeshT::facets_per_cell; ++j)
        x.insert( star[i].facet(mp, j).id(mp) );
    }
    // # dofs
    for (int i = 0; i < mapper.numVars(); ++i)
    {
      nnz += x.size()*mapper.variable(i).numDofsInFacet();
    }
  }

  // cells
  {
    // # dofs
    for (int i = 0; i < mapper.numVars(); ++i)
    {
      nnz += star.size()*mapper.variable(i).numDofsInCell();
    }
  }
  
  return nnz;
}

PetscErrorCode AppCtx::initPetscObjs()
{
//  printf("\nallocing petsc objs ... ");
  PetscErrorCode      ierr;

  for (unsigned i = 0; i < systems.size(); ++i)
  {
    System& sys = systems[i];

    for (unsigned j = 0; j < sys.Vec_u.size(); ++j)
    {
      ierr = VecCreate(PETSC_COMM_WORLD, &sys.Vec_u[j]);                     CHKERRQ(ierr);
      ierr = VecSetSizes(sys.Vec_u[j], PETSC_DECIDE, sys.n_dofs);            CHKERRQ(ierr);
      ierr = VecSetFromOptions(sys.Vec_u[j]);                                CHKERRQ(ierr);      
    }
    
    if (sys.active)
    {
      // create residue
      ierr = VecCreate(PETSC_COMM_WORLD, &sys.Vec_res);                     CHKERRQ(ierr);
      ierr = VecSetSizes(sys.Vec_res, PETSC_DECIDE, sys.n_dofs);            CHKERRQ(ierr);
      ierr = VecSetFromOptions(sys.Vec_res);                                CHKERRQ(ierr);
      
      // create jacobian
      ierr = MatCreate(PETSC_COMM_WORLD, &sys.Mat_jac);                                      CHKERRQ(ierr);
      ierr = MatSetSizes(sys.Mat_jac, PETSC_DECIDE, PETSC_DECIDE, sys.n_dofs, sys.n_dofs);   CHKERRQ(ierr);
      ierr = MatSetFromOptions(sys.Mat_jac);                                                 CHKERRQ(ierr);
      int nnz = estimateNnz(mp, *sys.mapper);
      //printf("NNZ = %d\n", nnz);
      //ierr = MatSeqAIJSetPreallocation(sys.Mat_jac, 0, nnz.data());                      CHKERRQ(ierr);
      ierr = MatSeqAIJSetPreallocation(sys.Mat_jac, nnz, PETSC_NULL);                      CHKERRQ(ierr);
      //ierr = MatSetOption(Mat_jac,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);                  CHKERRQ(ierr);
      ierr = MatSetOption(sys.Mat_jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);       CHKERRQ(ierr);
      
      // create nonlinear solver
      ierr = SNESCreate(PETSC_COMM_WORLD, &sys.snes);                                      CHKERRQ(ierr);
      ierr = SNESSetFunction(sys.snes, sys.Vec_res, FormFunction, this);                   CHKERRQ(ierr);
      ierr = SNESSetJacobian(sys.snes, sys.Mat_jac, sys.Mat_jac, FormJacobian, this);      CHKERRQ(ierr);
      //ierr = SNESSetJacobian(snes,Mat_jac,Mat_jac,SNESDefaultComputeJacobian,&user);  CHKERRQ(ierr);
     
      //ierr = SNESSetConvergenceTest(snes,CheckSnesConvergence,this,PETSC_NULL); CHKERRQ(ierr);
     
      ierr = SNESGetKSP(sys.snes,&sys.ksp);                                                  CHKERRQ(ierr);
      ierr = KSPGetPC(sys.ksp,&sys.pc);                                                      CHKERRQ(ierr);
      ierr = KSPSetOperators(sys.ksp,sys.Mat_jac,sys.Mat_jac,SAME_NONZERO_PATTERN);          CHKERRQ(ierr);
      //~ ierr = KSPSetType(ksp,KSPPREONLY);                                           CHKERRQ(ierr);
      //~ ierr = KSPSetType(ksp,KSPGMRES);                                               CHKERRQ(ierr);
      //~ ierr = PCSetType(pc,PCLU);                                                     CHKERRQ(ierr);
      //~ ierr = PCFactorSetMatOrderingType(pc, MATORDERINGNATURAL);                         CHKERRQ(ierr);
      //~ ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);  CHKERRQ(ierr);
    //~ ierr = SNESSetApplicationContext(snes,this);
     
      //~ #ifdef PETSC_HAVE_MUMPS
      //~ PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
      //~ #endif
     
      //ierr = SNESMonitorSet(snes, SNESMonitorDefault, 0, 0); CHKERRQ(ierr);
      //ierr = SNESMonitorSet(snes,Monitor,0,0);CHKERRQ(ierr);
      //ierr = SNESSetTolerances(snes,0,0,0,13,PETSC_DEFAULT);
//      ierr = SNESSetFromOptions(sys.snes); CHKERRQ(ierr);

      
    }
    
    
    
  }

  PetscFunctionReturn(0);
  //printf( "done\n");
}


#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
PetscErrorCode FormJacobian(SNES snes, Vec Vec_u, Mat *Mat_jac, Mat *prejac, MatStructure *flag, void *ptr)
{
//  AppCtx *user    = static_cast<AppCtx*>(ptr);
//  user->formJacobian(snes,Vec_u,Mat_jac,prejac,flag);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
PetscErrorCode FormFunction(SNES snes, Vec Vec_u, Vec Vec_fun, void *ptr)
{
//  AppCtx *user    = static_cast<AppCtx*>(ptr);
//  user->formFunction(snes,Vec_u,Vec_fun);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CheckSnesConvergence"
PetscErrorCode CheckSnesConvergence(SNES snes, PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
{
//  AppCtx *user    = static_cast<AppCtx*>(ctx);
//  user->checkSnesConvergence(snes, it, xnorm, pnorm, fnorm, reason);
  PetscFunctionReturn(0);
}
 





