#include "common.hpp"

extern void defineSystems(AppCtx& ctx);

void AppCtx::initAll()
{
  initShapes();
  defineSystems(*this);
  initMeshIo();
  initPetscObjs();
}

void AppCtx::initMeshIo()
{
  mesh_printer.attachMesh(mp);
}

void AppCtx::initShapes()
{
  shapes_c.resize(input.unk_fields.size());
  shapes_f.resize(input.unk_fields.size());
  shapes_r.resize(input.unk_fields.size());

  // defining the shape functions
  for (unsigned i = 0; i < shapes_c.size(); ++i)
  {
    shapes_c[i].setType( input.unk_fields[i].interpolation.c_str(), cell_dim(CELL_TYPE),   input.unk_fields[i].degree);
    shapes_f[i].setType( input.unk_fields[i].interpolation.c_str(), cell_dim(FACET_TYPE),  input.unk_fields[i].degree);
    if (cell_dim(CELL_TYPE) > 2)
      shapes_r[i].setType( input.unk_fields[i].interpolation.c_str(), cell_dim(RIDGE_TYPE),  input.unk_fields[i].degree);
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

PetscErrorCode AppCtx::initPetscObjs()
{
//  printf("\nallocing petsc objs ... ");
//  PetscErrorCode      ierr;
//
//  std::vector<UnkField> const& unk_fields = input.unk_fields;
//
//  Mat_jac.resize();
//  Vecs_u.resize();
//  Vecs_r.resize();
//  sness.resize();
//  ksps.resize();
//  pcs.resize();
//
//  for (unsigned sys = 0; sys < systems.size(); ++sys)
//  {
//    ierr = VecCreate(PETSC_COMM_WORLD, &Vec_res);                     CHKERRQ(ierr);
//    ierr = VecSetSizes(Vec_res, PETSC_DECIDE, n_unknowns);            CHKERRQ(ierr);
//    ierr = VecSetFromOptions(Vec_res);                                CHKERRQ(ierr);
//  }




  printf( "done\n");
}







