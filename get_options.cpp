#include "get_options.hpp"
#include <iostream>
#include <sstream>

std::string itos(int i)
{
  std::stringstream ss;
  ss << i;
  return ss.str();
}


// lua_states out in the same state
std::vector<double> callField(UserInput &user_in, int func_idx, double x, double y, double z, double t, int tag, std::vector<double> const& pars)
{
  lua_State *L = user_in.lua_state;
  int n_results = lua_gettop(L);

  lua_rawgeti(L, LUA_REGISTRYINDEX, func_idx);
  lua_pushnumber(L, x);
  lua_pushnumber(L, y);
  lua_pushnumber(L, z);
  lua_pushnumber(L, t);
  lua_pushnumber(L, tag);

  lua_newtable(L);
  for (int i = 0; i < (int)pars.size(); ++i)
  {
    lua_pushnumber(L, i); // push the index
    lua_pushnumber(L, pars[i]); // push the value at 'i'
    // the table is now at -3 on the stack. This tells Lua
    // to put "stack[-1]" at "stack[-2]" and pop them off,
    // leaving the table at the top of the stack
    lua_settable(L, -3);
  }

  /* do the call (6 arguments; LUA_MULTRET means unknown num of results) */
  if (lua_pcall(L, 6, LUA_MULTRET, 0) != 0)
  {
    lua_pushstring(L, "Attempt to call function(x,y,z,t,tag,pars).");
    lua_error(L);
    lua_close(L);
    throw;
  }

  n_results = lua_gettop(L) - n_results;

  if (n_results <= 0)
  {
    lua_pushstring(L, "Invalid number of results of function(x,y,z,t,tag,pars).");
    lua_error(L);
    lua_close(L);
    throw;
  }

  // retrieve results
  std::vector<double> results(n_results);
  for (int i = n_results-1; i >=0 ; --i)
  {
    results[i] =  luaL_checknumber(L,-1) ;
    lua_pop(L,1);
  }

  return results;
}

double callScalarField(UserInput &user_in, int func_idx, double x, double y, double z, double t, int tag, std::vector<double> const& pars)
{
  lua_State *L = user_in.lua_state;
  lua_rawgeti(L, LUA_REGISTRYINDEX, func_idx);
  lua_pushnumber(L, x);
  lua_pushnumber(L, y);
  lua_pushnumber(L, z);
  lua_pushnumber(L, t);
  lua_pushnumber(L, tag);

  lua_newtable(L);
  for (int i = 0; i < (int)pars.size(); ++i)
  {
    lua_pushnumber(L, i); // push the index
    lua_pushnumber(L, pars[i]); // push the value at 'i'
    // the table is now at -3 on the stack. This tells Lua
    // to put "stack[-1]" at "stack[-2]" and pop them off,
    // leaving the table at the top of the stack
    lua_settable(L, -3);
  }

  /* do the call (6 arguments; 1 result) */
  if (lua_pcall(L, 6, 1, 0) != 0)
  {
    lua_pushstring(L, "Attempt to call function(x,y,z,t,tag,pars).");
    lua_error(L);
    lua_close(L);
    throw;
  }

  double result =  luaL_checknumber(L,-1) ;
  lua_pop(L,1);

  return result;
}

void print_arg(int argc, const char* const* argv)
{
  int k = 0;
  for (int i = 0; i < argc; ++i)
  {
    while (argv[i][k] != '\0')
    {
      std::cout << argv[i][k];
      ++k;
    }
    std::cout << std::endl;
    k = 0;
  }

}

// throw if the top of the stack is not the correct type
#define MY_LUA_CHECK(type) void checkType_##type(lua_State *L, std::string const& name)            \
                           {                                                                       \
                             if (!lua_is##type(L, -1))                                             \
                             {                                                                     \
                               std::cout << "ERROR: LUA user input:\n";                            \
                               std::cout << "`" << name << "` should be a " << #type << std::endl; \
                               lua_pushstring(L, "incorrect argument");                            \
                               lua_error(L);                                                       \
                               lua_close(L);                                                       \
                               throw;                                                              \
                             }                                                                     \
                           }

MY_LUA_CHECK(string)
MY_LUA_CHECK(table)
MY_LUA_CHECK(number)
MY_LUA_CHECK(boolean)
MY_LUA_CHECK(none)
MY_LUA_CHECK(function)
MY_LUA_CHECK(cfunction)
MY_LUA_CHECK(lightuserdata)
MY_LUA_CHECK(userdata)

#undef MY_LUA_CHECK

bool read_user_settings(lua_State *lua_state, UserInput::Settings & settings)
{
  lua_getglobal(lua_state, "settings");
  checkType_table(lua_state, "settings");


  lua_getfield(lua_state, -1, "name");
  checkType_string(lua_state, "name");
  settings.name = lua_tostring(lua_state, -1);
  lua_pop(lua_state, 1);  // pops name

  lua_getfield(lua_state, -1, "space_dim");
  checkType_number(lua_state, "space_dim");
  settings.space_dim = lua_tonumber(lua_state, -1);
  lua_pop(lua_state, 1);  // pops space_dim


  lua_getfield(lua_state, -1, "times");
  checkType_table(lua_state, "times");
  {
    lua_getfield(lua_state, -1, "from");
    checkType_number(lua_state, "from");
    settings.time_from = lua_tonumber(lua_state, -1);
    lua_pop(lua_state, 1); // pops from

    lua_getfield(lua_state, -1, "to");
    checkType_number(lua_state, "to");
    settings.time_to = lua_tonumber(lua_state, -1);
    lua_pop(lua_state, 1); // pops to

    lua_getfield(lua_state, -1, "step");
    checkType_number(lua_state, "step");
    settings.time_step = lua_tonumber(lua_state, -1);
    lua_pop(lua_state, 1); // pops step

    if (settings.time_to < settings.time_from)
    {
      std::cout << "ERROR LUA USER: final time less than initial time\n";
      lua_close(lua_state);
      throw;
    }
    if (settings.time_step < 0)
    {
      std::cout << "ERROR LUA USER: negative time step\n";
      lua_close(lua_state);
      throw;
    }
    settings.time_current = settings.time_from;
  }
  lua_pop(lua_state, 1);  // pops times


  lua_getfield(lua_state, -1, "stop_at_steady_state");
  checkType_boolean(lua_state, "stop_at_steady_state");
  settings.stop_at_steady_state = lua_toboolean(lua_state, -1);
  lua_pop(lua_state, 1);  // pops stop_at_steady_state

  if (settings.stop_at_steady_state)
  {
    lua_getfield(lua_state, -1, "steady_state_tolerance");
    checkType_number(lua_state, "steady_state_tolerance");
    settings.steady_tol = lua_tonumber(lua_state, -1);
    lua_pop(lua_state, 1);  // pops steady_state_tolerance
  }

  //std::cout << settings.name << std::endl;
  //std::cout << settings.space_dim << std::endl;

  lua_getfield(lua_state, -1, "print_step");
  checkType_number(lua_state, "print_step");
  settings.print_step = lua_tonumber(lua_state, -1);
  lua_pop(lua_state, 1); // pops print_step

  lua_getfield(lua_state, -1, "degree_of_pol_integrated_exactly_for_cells");
  checkType_number(lua_state, "degree_of_pol_integrated_exactly_for_cells");
  settings.quadr_degree_cell = lua_tonumber(lua_state, -1);
  lua_pop(lua_state, 1); // pops degree_of_pol_integrated_exactly_for_cells

  lua_getfield(lua_state, -1, "degree_of_pol_integrated_exactly_for_facets");
  checkType_number(lua_state, "degree_of_pol_integrated_exactly_for_facets");
  settings.quadr_degree_facet = lua_tonumber(lua_state, -1);
  lua_pop(lua_state, 1); // pops degree_of_pol_integrated_exactly_for_facets

  lua_getfield(lua_state, -1, "degree_of_pol_integrated_exactly_for_ridges");
  checkType_number(lua_state, "degree_of_pol_integrated_exactly_for_ridges");
  settings.quadr_degree_ridge = lua_tonumber(lua_state, -1);
  lua_pop(lua_state, 1); // pops degree_of_pol_integrated_exactly_for_ridges

  lua_getfield(lua_state, -1, "degree_of_pol_integrated_exactly_for_errors");
  checkType_number(lua_state, "degree_of_pol_integrated_exactly_for_errors");
  settings.quadr_degree_err = lua_tonumber(lua_state, -1);
  lua_pop(lua_state, 1); // pops degree_of_pol_integrated_exactly_for_errors

  lua_getfield(lua_state, -1, "print_errors");
  checkType_boolean(lua_state, "print_errors");
  settings.print_errors = lua_toboolean(lua_state, -1);
  lua_pop(lua_state, 1); // pops `print_errors`

  lua_getfield(lua_state, -1, "integ_quantities_names");
  checkType_table(lua_state, "integ_quantities_names");
  for (int name_i = 1; ; ++name_i)
  {
    lua_rawgeti(lua_state, -1, name_i);
    if (lua_isnil(lua_state,-1))
    {
      lua_pop(lua_state,1);
      break;
    }

    checkType_string(lua_state, "integ_quantities_names[i]");
    settings.integ_quant_names.push_back(lua_tostring(lua_state, -1));

    lua_pop(lua_state, 1); // pops `name_i`
  }
  lua_pop(lua_state, 1); // pops `integ_quantities_names`

  lua_getfield(lua_state, -1, "print_integ_quantities");
  checkType_table(lua_state, "print_integ_quantities");
  for (int p_i = 1; ; ++p_i)
  {
    lua_rawgeti(lua_state, -1, p_i);
    if (lua_isnil(lua_state,-1))
    {
      lua_pop(lua_state,1);
      break;
    }

    checkType_boolean(lua_state, "print_integ_quantities[i]");
    settings.integ_quant_print.push_back(lua_toboolean(lua_state, -1));

    lua_pop(lua_state, 1); // pops `p_i`
  }
  lua_pop(lua_state, 1); // pops `print_integ_quantities`

  lua_getfield(lua_state, -1, "input_meshfile");
  checkType_string(lua_state, "input_meshfile");
  settings.input_meshfile = lua_tostring(lua_state, -1);
  lua_pop(lua_state, 1); // pops input_meshfile

  lua_getfield(lua_state, -1, "output_basename");
  checkType_string(lua_state, "output_basename");
  settings.output_basename = lua_tostring(lua_state, -1);
  lua_pop(lua_state, 1); // pops output_basename

  lua_pop(lua_state, 1); // pop settings

  return 0;
}

bool read_user_region_marks(lua_State *lua_state, UserInput::RegionMarks & regions)
{
  lua_getglobal(lua_state, "region_marks");
  checkType_table(lua_state, "region_marks");

  //
  // Reading defined regions
  //
  for (int reg = 1; ; ++reg)
  {
    lua_rawgeti(lua_state, -1, reg);
    if (lua_isnil(lua_state,-1))
    {
      lua_pop(lua_state,1);
      break;
    }

    lua_getfield(lua_state, -1, "name");
    checkType_string(lua_state, "name");
    regions.names.push_back( lua_tostring(lua_state, -1) );
    lua_pop(lua_state, 1);

    regions.tags.push_back(std::vector<int>());
    lua_getfield(lua_state, -1, "tags");
    {
      checkType_table(lua_state, "tags");
      for (int jj = 1; ; ++jj)
      {
        lua_rawgeti(lua_state, -1, jj);
        if (lua_isnil(lua_state,-1))
        {
          lua_pop(lua_state,1);
          break;
        }
        checkType_number(lua_state, "tags[]");

        regions.tags.back().push_back(lua_tonumber(lua_state,-1));

        lua_pop(lua_state,1); // pop dirich_tags[i][j]
      }
    }
    lua_pop(lua_state, 1); // pop tags

    lua_pop(lua_state, 1); // pop reg
  }

  lua_getfield(lua_state, -1, "do_after_reading");
  checkType_string(lua_state, "do_after_reading");
  regions.do_after_reading = lua_tostring(lua_state,-1);
  //std::cout << regions.do_after_reading << std::endl;
  lua_pop(lua_state, 1); // pops do_after_reading


  lua_pop(lua_state, 1); // pops region_marks

  return 0;
}

bool read_user_unknown_fields(lua_State *lua_state, std::vector<UserInput::UnkField> & unk_fields)
{
  lua_getglobal(lua_state, "unknown_fields");
  checkType_table(lua_state, "unknown_fields");

  unk_fields.clear();

  // looping over fields
  for (int i_unk = 1; ; ++i_unk)
  {
    lua_rawgeti(lua_state, -1, i_unk);
    if (lua_isnil(lua_state,-1))
    {
      lua_pop(lua_state,1);
      break;
    }

    unk_fields.push_back(UserInput::UnkField());

    UserInput::UnkField & unk_i = unk_fields.back();

    // for each field, read its properties
    checkType_table(lua_state, "unknown_fields[i_unk]");
    {
      lua_getfield(lua_state, -1, "name");
      checkType_string(lua_state, "name");
      unk_i.name = lua_tostring(lua_state, -1);
      lua_pop(lua_state, 1);

      lua_getfield(lua_state, -1, "n_components");
      checkType_number(lua_state, "n_components");
      unk_i.n_comps = lua_tonumber(lua_state, -1);
      lua_pop(lua_state, 1);

      lua_getfield(lua_state, -1, "degree");
      checkType_number(lua_state, "degree");
      unk_i.degree = lua_tonumber(lua_state, -1);
      lua_pop(lua_state, 1);

      lua_getfield(lua_state, -1, "interpolation");
      checkType_string(lua_state, "interpolation");
      unk_i.interpolation = lua_tostring(lua_state, -1);
      lua_pop(lua_state, 1);

      // Where is defined
      //
      lua_getfield(lua_state, -1, "where_is_defined");
      {
        checkType_table(lua_state, "where_is_defined");

        // loop on components
        for (int cp = 1; ; ++cp)
        {
          lua_rawgeti(lua_state, -1, cp);
          if (lua_isnil(lua_state,-1))
          {
            lua_pop(lua_state,1);
            break;
          }
          checkType_table(lua_state, "where_is_defined");

          unk_i.where_is_defined.push_back(std::vector<int>());

          for (int jj = 1; ; ++jj)
          {
            lua_rawgeti(lua_state, -1, jj);
            if (lua_isnil(lua_state,-1))
            {
              lua_pop(lua_state,1);
              break;
            }
            checkType_number(lua_state, "where_is_defined[]");

            unk_i.where_is_defined.back().push_back(lua_tonumber(lua_state,-1));

            lua_pop(lua_state,1); // pop where_is_defined[i][j]
          }


          lua_pop(lua_state,1);  // pop where_is_defined[i]
        }

      }
      lua_pop(lua_state, 1); // pop where_is_defined
      // decompose where_is_defined
      unk_i.n_regions = unk_i.where_is_defined.size();
      unk_i.n_tags.resize(unk_i.n_regions);
      unk_i.tags.clear();
      for (int k = 0; k < unk_i.n_regions; ++k)
      {
        unk_i.n_tags[k] = unk_i.where_is_defined[k].size();
        for (int j = 0; j < unk_i.n_tags[k]; ++j)
          unk_i.tags.push_back( unk_i.where_is_defined[k][j] );
      }



      // Dirichlet tags
      //
      lua_getfield(lua_state, -1, "dirichlet_tags");
      {
        checkType_table(lua_state, "dirichlet_tags");

        for (int jj = 1; ; ++jj)
        {
          lua_rawgeti(lua_state, -1, jj);
          if (lua_isnil(lua_state,-1))
          {
            lua_pop(lua_state,1);
            break;
          }
          checkType_number(lua_state, "dirichlet_tags[]");

          unk_i.dirich_tags.push_back(lua_tonumber(lua_state,-1));

          lua_pop(lua_state,1); // pop dirich_tags[i][j]
        }
      }
      lua_pop(lua_state, 1); // pop dirichlet_tags


      //lua_getfield(lua_state, -1, "exact_val");
      //checkType_function(lua_state, "exact_val");
      //unk_i.exact_value_idx = luaL_ref(lua_state, LUA_REGISTRYINDEX);
      // it doesn't need lua_pop .. lua_ref already pops

    }

    lua_pop(lua_state, 1); // pop i_unk
  } // end loop on unknowns


  lua_pop(lua_state, 1); // pop unknown_fields

  return 0;
}

bool read_user_systems(lua_State *lua_state, std::vector<UserInput::System>& systems)
{
  lua_getglobal(lua_state, "systems");
  checkType_table(lua_state, "systems");

  systems.clear();

  // looping over systems
  for (int sys_i = 1; ; ++sys_i)
  {
    lua_rawgeti(lua_state, -1, sys_i);
    if (lua_isnil(lua_state,-1))
    {
      lua_pop(lua_state,1);
      break;
    }

    systems.push_back(UserInput::System());
    UserInput::System& sys = systems.back();

    // for each system, read its properties
    checkType_table(lua_state, "systems[sys]");
    {
      lua_getfield(lua_state, -1, "actived");
      checkType_boolean(lua_state, "actived");
      sys.actived = lua_toboolean(lua_state, -1);
      lua_pop(lua_state, 1);

      lua_getfield(lua_state, -1, "first_dof_number");
      checkType_number(lua_state, "first_dof_number");
      sys.first_dof_number = lua_tonumber(lua_state, -1);
      lua_pop(lua_state, 1);

      lua_getfield(lua_state, -1, "fields");
      {
        checkType_table(lua_state, "fields");

        for (int jj = 1; ; ++jj)
        {
          lua_rawgeti(lua_state, -1, jj);
          if (lua_isnil(lua_state,-1))
          {
            lua_pop(lua_state,1);
            break;
          }
          checkType_number(lua_state, "fields[]");

          sys.fields_id.push_back(lua_tonumber(lua_state,-1));

          lua_pop(lua_state,1); // pop fields[i]
        }
      }
      lua_pop(lua_state, 1); // pop fields

      lua_getfield(lua_state, -1, "n_copies_per_ts");
      checkType_number(lua_state, "n_copies_per_ts");
      sys.n_copies_per_ts = lua_tonumber(lua_state, -1);
      lua_pop(lua_state, 1);

      lua_getfield(lua_state, -1, "mapper");
      if (lua_isnumber(lua_state, -1))
      {
        int a = lua_tonumber(lua_state, -1);
        if (a<0)
          { printf("\nERROR:invalid mapper\n"); throw;}
        sys.mapper = a;
      }
      lua_pop(lua_state, 1);

    }


    lua_pop(lua_state, 1); // pop sys_i
  }


  lua_pop(lua_state, 1); // pop systems

  return 0;
}

bool read_user_parameters(lua_State *lua_state, std::vector<UserInput::Parameter> & parameters)
{
  lua_getglobal(lua_state, "parameters");
  checkType_table(lua_state, "parameters");

  for (int idx = 1; ; ++idx)
  {
    lua_rawgeti(lua_state, -1, idx);

    if (lua_isnil(lua_state,-1))
    {
      lua_pop(lua_state,1);
      break;
    }

    parameters.push_back(UserInput::Parameter());
    //parameters.back().name = lua_tostring(lua_state,-2);

    checkType_table(lua_state, std::string("parameters[")+itos(idx-1)+std::string("]"));
    {
      lua_getfield(lua_state, -1, "name");
      checkType_string(lua_state, "name");
      parameters.back().name = lua_tostring(lua_state,-1);
      lua_pop(lua_state,1);

      lua_getfield(lua_state, -1, "values");
      int values_index = lua_gettop(lua_state);
      if (lua_isnumber(lua_state,-1))
        parameters.back().values.push_back(lua_tonumber(lua_state,-1));
      else
      if (lua_istable(lua_state,-1))
      {
        lua_getfield(lua_state, -1, "from");
        if (lua_isnumber(lua_state,-1))  // arithmetic progression values
        {
          double from = 0, to = 0, step = 0;
          from = lua_tonumber(lua_state,-1);
          lua_pop(lua_state,1);

          lua_getfield(lua_state, -1, "to");
          checkType_number(lua_state, "to");
          to = lua_tonumber(lua_state,-1);
          lua_pop(lua_state,1);

          lua_getfield(lua_state, -1, "step");
          checkType_number(lua_state, "step");
          step = lua_tonumber(lua_state,-1);
          lua_pop(lua_state,1);

          if (to >= from && step > 0)
          {
            for (double a = from; a <= to; a += step)
              parameters.back().values.push_back(a);
          }
          else if (to <= from && step < 0)
          {
            for (double a = from; a >= to; a += step)
              parameters.back().values.push_back(a);
          }
          else {
            std::cout << "ERROR LUA USER: infinite parameters values\n";
            lua_close(lua_state);
            throw;
          }
        }
        else
        {
          lua_pop(lua_state,-1);
          lua_settop(lua_state, values_index);
          for (int kk = 1; ;++kk)
          {
            lua_rawgeti(lua_state, -1, kk);
            if (lua_isnil(lua_state,-1))
            {
              lua_pop(lua_state,1);
              break;
            }
            checkType_number(lua_state, "values[]");
            parameters.back().values.push_back(lua_tonumber(lua_state,-1));
            lua_pop(lua_state,1);
          }
        }
      }
      else {
        std::cout << "ERROR LUA USER: invalid parameter value\n";
        lua_close(lua_state);
        throw;
      }
      lua_pop(lua_state,1);
    }

    lua_pop(lua_state, 1);
  }

  lua_pop(lua_state, 1); // pop parameters

  return 0;
}

bool read_user_petsc_options(lua_State *lua_state, UserInput::LuaPetscOps & petsc_opts)
{

  lua_getglobal(lua_state, "petsc_options");
  checkType_table(lua_state, "petsc_options");

  lua_rawgeti(lua_state, -1, 1);

  checkType_string(lua_state, "petsc_options[1]");

  std::string opts(luaL_checkstring(lua_state, -1));
  opts = "MAIN_PROGRAM\n" + opts;

  std::stringstream iss(opts);
  std::string token;
  while( iss >> token)
  {
    char *arg = new char[token.size() + 1];
    for (unsigned i = 0; i < token.size(); ++i)
      arg[i] = token[i];
    arg[token.size()] = '\0';
    petsc_opts.tokens.push_back(arg);
  }

  // debug
  //print_arg(petsc_opts.tokens.size(), petsc_opts.tokens.data());

  lua_pop(lua_state, 1); // pops petsc_options[1]

  lua_pop(lua_state, 1); // pops petsc_options

  return 0;
}

int read_input(UserInput & user_in)
{
  //std::vector<System>    &  systems      = user_in.systems     ;
  std::vector<UserInput::UnkField>  &  unk_fields   = user_in.unk_fields  ;
  std::vector<UserInput::System>    &  systems      = user_in.systems     ;
  UserInput::RegionMarks            &  regions      = user_in.regions     ;
  UserInput::Settings               &  settings     = user_in.settings    ;
  std::vector<UserInput::Parameter> &  parameters   = user_in.parameters  ;
  UserInput::LuaPetscOps            &  petsc_opts   = user_in.petsc_opts  ;
  lua_State*             &  lua_state    = user_in.lua_state   ;
  std::string            &  cfg_fname    = user_in.cfg_fname   ;

  lua_state = luaL_newstate();

  // Opens all standard Lua libraries into the given state.
  luaL_openlibs(lua_state);


  // read file
  if ( luaL_dofile(lua_state, cfg_fname.c_str()) )
  {
    std::cout << "ERROR: invalid file format" << std::endl;
    std::cout << lua_tostring(lua_state, -1) << std::endl;

    lua_close(lua_state);
    return -1;
  }



  // --------------------------------------------------------------
  // UNKNOWN FIELDS
  // --------------------------------------------------------------

  if (read_user_unknown_fields(lua_state, unk_fields))
    throw;

  // --------------------------------------------------------------
  // SYSTEMS FIELDS
  // --------------------------------------------------------------

  if (read_user_systems(lua_state, systems))
    throw;

  // --------------------------------------------------------------
  // REGION MARKS
  // --------------------------------------------------------------

  if (read_user_region_marks(lua_state, regions))
    throw;


  // --------------------------------------------------------------
  // SETTINGS
  // --------------------------------------------------------------

  if (read_user_settings(lua_state, settings) )
    throw;

  // --------------------------------------------------------------
  // PARAMETERS
  // --------------------------------------------------------------

  if(read_user_parameters(lua_state, parameters))
    throw;

  //
  //  PETSC_OPTIONS
  //

  if (read_user_petsc_options(lua_state, petsc_opts))
    throw;

  return 0;
}





