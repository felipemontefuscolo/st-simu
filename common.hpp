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
#include <lua.hpp>
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


// lua_states out in the same state
std::vector<double> callField(AppCtx& user_in, int func_idx, double x, double y, double z, double t, int tag, std::vector<double> const& pars);
double callScalarField(AppCtx& user_in, int func_idx, double x, double y, double z, double t, int tag, std::vector<double> const& pars);
int read_input(AppCtx& user_in);



struct AppCtx
{
  struct Settings
  {


    std::string name;
    int space_dim;

    double time_from;
    double time_to;
    double time_step;

    bool stop_at_steady_state;
    double steady_tol;

    int                      print_step;
    int                      quadr_degree_cell;
    int                      quadr_degree_facet;
    int                      quadr_degree_ridge;
    int                      quadr_degree_err;
    bool                     print_errors;
    std::vector<std::string> integ_quant_names;
    std::vector<bool>        integ_quant_print;
    std::vector<double>      integ_quant_vals;
    std::string              input_meshfile;
    std::string              output_basename;
  };

  struct UnkField
  {
    UnkField() : name("undefined"), n_comps(1), sys_num(-1) {}
    
    std::string name;
    int n_comps;
    std::string interpolation;
    int degree;
    int sys_num;

    std::vector<std::vector<int> > where_is_defined; // each where_is_defined[i] is a set of tags
    // decomposition of the where_is_defined vector
    int n_regions;
    std::vector<int> n_tags;
    std::vector<int> tags;


    std::vector<int> dirich_tags; //

    //int exact_value_idx; // index in the LUA registry where the functions is.
  };

  struct System
  {
    System() : fields_id(), first_dof_number(0), n_copies_per_ts(2), sys_mapper(-1), active(true) {}

    std::vector<int> fields_id; // unk_fields id of this system
    int first_dof_number;
    int n_copies_per_ts;
    int sys_mapper; // the id of another system from which the mapper is copied; -1 means create its own mapper;
    bool active;

    // Petsc structures
    std::vector<Vec> Vec_u; // Vec_u.size() = n_copies_per_ts; Vec_u[0] is the u at time n+1, Vec_u[1] at n, and so on
    Vec Vec_res; // it will be created only if active=true
    Mat Mat_jac; // it will be created only if active=true
    SNES snes;   // it will be created only if active=true
    KSP ksp;     // it will be created only if active=true
    PC pc;       // it will be created only if active=true

    // Mapper: maps the mesh to the degrees of freedom
    std::tr1::shared_ptr<DofMapperT> mapper; // if sys_mapper < 0, it's allocated a new mapper; otherwise this->mapper=sys_mapper.mapper
    alelib::index_t n_dofs; // warning: In a more general case this is different from mapper.numDofs()
  };

  struct RegionMarks
  {
    //RegionMarks() : solid_tags(),gas_tags(),cl_tags(),do_after_reading() {}

    std::vector< std::string >      names;
    std::vector< std::vector<int> > tags;

    std::string do_after_reading;
  };

  struct Parameter
  {
    Parameter() : current_position(0)
    { }
    std::string name;
    std::vector<double> values;
    size_t current_position;
    double val() const
    { return values.at(current_position); }
  };

  struct LuaPetscOps
  {
    LuaPetscOps() {}

    std::vector<char* > tokens; // options splitted in tokens: same as argv

    ~LuaPetscOps()
    {
      for (unsigned i = 0; i < tokens.size(); ++i)
      {
        if (tokens[i])
          delete [] tokens[i];
      }
    }
  };

  AppCtx() : lua_state(0), mp(NULL)
  { }

  inline ~AppCtx(); // defined below

  //
  // User input members
  //
  Settings                 settings;
  std::vector<UnkField>    unk_fields;
  std::vector<System>      systems;
  RegionMarks              regions;

  std::vector<Parameter>   parameters;
  LuaPetscOps              petsc_opts;
  lua_State                *lua_state;
  std::string              cfg_fname;
  //
  // End User input
  //

  Timer timer;

  MeshT* mp;
  
  int                 current_system_to_solve;
  double              current_time;
  std::vector<double> current_pars; // parameters

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






  // ******************** MEMBER FUNCTIONS *************************

  void printSettings() const;


  void initAll();
  void initShapes();
  void initMesh();
  void initMeshIo();
  void initDofMappers();
  PetscErrorCode initPetscObjs();
  PetscErrorCode formJacobian(SNES /*snes*/,Vec x,Mat *Mat_Jac, Mat* /*prejac*/, MatStructure * /*flag*/);
  PetscErrorCode formFunction(SNES /*snes*/, Vec x, Vec f);

  // solving the systems
  void solveTemporalProblem();
  void setUpInitialConditions();




  void destroyPetsc();
};


AppCtx::~AppCtx()
{
  if (mp)
    delete mp;
  /* cleanup Lua */
  if (lua_state != NULL)
  {
    int n = lua_gettop(lua_state);
    printf("\nNumber of elements in lua stack : %i\n", n);
    lua_close(lua_state);
  }

}







#endif




