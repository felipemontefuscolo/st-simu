#ifndef GET_OPTIONS_HPP
#define GET_OPTIONS_HPP

#include <vector>
#include <string>
#include <cstdio>
#include <lua.hpp>


struct UserInput;

// lua_states out in the same state
std::vector<double> callField(UserInput &user_in, int func_idx, double x, double y, double z, double t, int tag, std::vector<double> const& pars);


double callScalarField(UserInput &user_in, int func_idx, double x, double y, double z, double t, int tag, std::vector<double> const& pars);




struct UserInput
{
  struct Settings
  {
    
    
    std::string name;
    int space_dim;

    double time_from;
    double time_to;
    double time_step;
    double time_current;

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
    System() : fields_id(), first_dof_number(0), n_copies_per_ts(2), mapper(-1), actived(true) {}
    
    std::vector<int> fields_id; // unk_fields id of this system
    int first_dof_number;
    int n_copies_per_ts;
    int mapper; // the id of another system from which the mapper is copied; -1 means create its own mapper;
    bool actived;
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
    
  
  UserInput() : lua_state(0) {}

  Settings                 settings;
  std::vector<UnkField>    unk_fields;
  std::vector<System>      systems;
  RegionMarks              regions;

  std::vector<Parameter>   parameters;
  LuaPetscOps              petsc_opts;
  lua_State                *lua_state;
  std::string              cfg_fname;
  
  ~UserInput()
  {
    /* cleanup Lua */
    if (lua_state != NULL)
    {
      int n = lua_gettop(lua_state);
      printf("\nNumber of elements in lua stack : %i\n", n);
      lua_close(lua_state);
    }
  }
};


int read_input(UserInput & user_in);

#endif



