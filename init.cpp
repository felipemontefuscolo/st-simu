#include "common.hpp"


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

  cout << endl;
  cout << "Mesh data:" << endl;  
  cout << "==========" << endl;  
  cout << "# vertices   : " << mp->numVertices() << endl;
  cout << "# ridges(3D) : " << mp->numRidges() << endl;
  cout << "# facets     : " << mp->numFacets() << endl;
  cout << "# cells      : " << mp->numCells() << endl;
  
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

  cout << endl;
  cout << "Degrees of freedom:" << endl;
  cout << "===================" << endl;
  cout << "# dofs:" << endl;
  for (unsigned i = 0; i < systems.size(); ++i)
    if (systems[i].active)
      cout << "(active) system (" << i << ") : " << systems[i].n_dofs << endl;
    else
      cout << "         system (" << i << ") : " << systems[i].n_dofs << endl;
  
}

PetscErrorCode AppCtx::initPetscObjs()
{
//  printf("\nallocing petsc objs ... ");
  PetscErrorCode      ierr;

  for (unsigned i = 0; i < systems.size(); ++i)
  {
    System& sys = systems[i];
    ierr = VecCreate(PETSC_COMM_WORLD, &sys.Vec_res);                     CHKERRQ(ierr);
    ierr = VecSetSizes(sys.Vec_res, PETSC_DECIDE, sys.n_dofs);            CHKERRQ(ierr);
    ierr = VecSetFromOptions(sys.Vec_res);                                CHKERRQ(ierr);
    
    //for (unsigned j = 0; j < sys.Vec_u.size(); ++j)
    //{
    //  ierr = VecCreate(PETSC_COMM_WORLD, &sys.Vec_u[j]);                     CHKERRQ(ierr);
    //  ierr = VecSetSizes(sys.Vec_u[j], PETSC_DECIDE, sys.n_dofs);            CHKERRQ(ierr);
    //  ierr = VecSetFromOptions(sys.Vec_u[j]);                                CHKERRQ(ierr);      
    //}
  }

  PetscFunctionReturn(0);
  //printf( "done\n");
}







