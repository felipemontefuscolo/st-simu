#ifndef COMMMON_HPP
#define COMMMON_HPP


#include <Alelib/Mesh>
#include <Alelib/Quadrature>
#include <Alelib/DofMapper>
#include <Alelib/ShapeFunction>
#include <Alelib/IO>
#include "petscsnes.h"
#include <iostream>
#include <tr1/memory>
#include <tr1/array>
#include "get_options.hpp"
//#include "mypetsc.hpp"

#if   CELL_TYPE==TRIANGLE
# define FACET_TYPE EDGE
# define RIDGE_TYPE POINT
#elif CELL_TYPE==TETRAHEDRON
# define FACET_TYPE TRIANGLE
# define RIDGE_TYPE EDGE
#elif CELL_TYPE==EDGE
# define FACET_TYPE POINT
# define RIDGE_TYPE POINT
#elif CELL_TYPE==QUADRANGLE
# define FACET_TYPE EDGE
# define RIDGE_TYPE POINT
#elif CELL_TYPE==HEXAHEDRON
# define FACET_TYPE QUADRANGLE
# define RIDGE_TYPE EDGE
#endif

using namespace std;
using namespace alelib;
using namespace tr1;

// define CELL_TYPE in the Makefile

typedef Mesh<CELL_TYPE> MeshT;
typedef typename MeshT::VertexH VertexH;
typedef typename MeshT::RidgeH  RidgeH;
typedef typename MeshT::FacetH  FacetH;
typedef typename MeshT::CellH   CellH;

typedef DofMapper<MeshT>       DofMapperT;
typedef MeshIoMsh<CELL_TYPE>   IoMshT;
typedef MeshIoVtk<CELL_TYPE>   IoVtkT;

typedef marray::Array<double, 1, marray::RowMajor, double[3]> Vector;
typedef marray::Array<double, 2, marray::RowMajor, double[9]> Tensor;


const double pi  = 3.141592653589793;
const double pi2 = pi*pi;
inline double sqr(double v) {return v*v;}

template<class Vec, class T>
bool is_in(T value, Vec const& v)
{
  for (int i = 0; i < (int)v.size(); i++)
  {
    if (value==v[i])
      return true;
  }
  return false;
}

struct AppCtx;


// tricks to avoid compiler error about OPENMP
#ifndef ALE_HAS_OPENMP

inline int omp_get_thread_num() {return 0;};
inline int omp_get_num_threads(){return 1;};

#endif

static inline
int cell_dim(ECellType ct)
{
  switch (ct)
  {
    case TRIANGLE:
      return 2;
    case TETRAHEDRON:
      return 3;
    case EDGE:
      return 1;
    case QUADRANGLE:
      return 2;
    case HEXAHEDRON:
      return 3;
    default: throw;
  }
}

struct AppCtx
{
  AppCtx() : mp(NULL) {}

  UserInput input;

  Timer timer;

  MeshT* mp;

  std::vector<ShapeFunction>   shapes_c; // cell ... shapes_c[i] returns the shape-function of the variable i
  std::vector<ShapeFunction>   shapes_f; // facet
  std::vector<ShapeFunction>   shapes_r; // ridge

  // shape functions evaluated at quadrature points
  marray::Array<std::vector<double>, 2> phi_c; // phi[qp][var][i-th]
  marray::Array<std::vector<double>, 2> phi_f;
  marray::Array<std::vector<double>, 2> phi_r;

  // derivative of the shape functions evaluated at quadrature points (master element)
  marray::Array<marray::Array<double,2>, 2> dLphi_c; // dphi[qp][var][i-th][comp]
  marray::Array<marray::Array<double,2>, 2> dLphi_f;
  marray::Array<marray::Array<double,2>, 2> dLphi_r;

  std::vector<DofMapperT> dof_mappers;

  Quadrature Q_c;
  Quadrature Q_f;
  Quadrature Q_r;

  IoMshT   mesh_reader;
  IoVtkT   mesh_printer;
  MeshT    mesh;

  // PETSC OBJS
  std::vector<Mat>  Mat_jac; // vector jacobians
  std::vector<Vec>  Vecs_u; // vector of unknowns vector
  std::vector<Vec>  Vecs_r; // vector of residues
  std::vector<SNES> sness; // vector of non-linear system context
  std::vector<KSP>  ksps;
  std::vector<PC>   pcs;
  

  void initAll();
  void initShapes();
  void initMeshIo();
  PetscErrorCode initPetscObjs();

  ~AppCtx()
  {
    if (mp)
      delete mp;
  }
};










#endif




