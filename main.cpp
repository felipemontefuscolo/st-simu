#include "common.hpp"

static char help[] = "\nUsage: %s <path_to_your_cfg_file.lua> \n\nExample:\n %s some_dir/test.lua\n\n";

const char* cell_name(ECellType ct)
{
  switch (ct)
  {
    case TRIANGLE:
      return "TRIANGLE";
    case TETRAHEDRON:
      return "TETRAHEDRON";
    case EDGE:
      return "EDGE";
    case QUADRANGLE:
      return "QUADRANGLE";
    case HEXAHEDRON:
      return "HEXAHEDRON";
    default: throw;
  }
}

void print_settings(AppCtx & user_ctx)
{
  //cout << "Cell type number " << CELL_TYPE << endl;
  //cout << "Cell type number " << FACET_TYPE << endl;
  //cout << "Cell type number " << RIDGE_TYPE << endl;
  
  cout << "Settings:\n";
  cout << "---------\n";
  cout << "name                        : " << user_ctx.input.settings.name << endl;
  cout << "mesh type                   : " << cell_name(CELL_TYPE) << endl;
  cout << "space dimension             : " << user_ctx.input.settings.space_dim << endl;
  cout << "time                        : " << "from=" << user_ctx.input.settings.time_from
                                           <<", to="<<user_ctx.input.settings.time_to
                                           <<", step=" << user_ctx.input.settings.time_step << endl;
  cout << "exact integ. degree (cell)  : " << user_ctx.input.settings.quadr_degree_cell << endl;
  cout << "exact integ. degree (facet) : " << user_ctx.input.settings.quadr_degree_facet << endl;
  cout << "exact integ. degree (ridge) : " << user_ctx.input.settings.quadr_degree_ridge << endl;
  cout << "input mesh file             : " << user_ctx.input.settings.input_meshfile << endl;
  cout << "output mesh file (basename) : " << user_ctx.input.settings.output_basename << endl;

  cout << endl;

  cout << "Variables:\n";
  cout << "----------\n";
  {
    for (int j = 0; j < (int)user_ctx.input.unk_fields.size(); ++j)
    {
      UserInput::UnkField & var = user_ctx.input.unk_fields[j];
      cout << "name              : " << var.name << endl;
      cout << "num of components : " << var.n_comps << endl;
      cout << "interpolation     : " << var.interpolation << endl;
      cout << "degree            : " << var.degree << endl;
      cout << "system #          : " << var.sys_num << endl;
      cout << "where is defined  : " << "{ ";
      for (unsigned i = 0; i < var.where_is_defined.size(); ++i)
      {
        cout << "{";
        for (unsigned k = 0; k < var.where_is_defined[i].size(); ++k)
        {
          cout << var.where_is_defined[i][k];
          if (k != var.where_is_defined[i].size() - 1)
            cout << ",";
        }
        cout << "}";
        if (i != var.where_is_defined.size() - 1)
          cout << ", ";
        else
          cout << " }" << endl;
      }

      cout << "dirichlet tags    : " << "{ ";
      for (unsigned i = 0; i < var.dirich_tags.size(); ++i)
      {
        cout << var.dirich_tags[i];
        if (i != var.dirich_tags.size()-1)
          cout << ",";
        else
          cout << " }" << endl;
      }
      cout << endl;
    }
    cout << endl;
  }

  cout << "Regions:\n";
  cout << "--------\n";
  for (unsigned i = 0; i < user_ctx.input.regions.names.size(); ++i)
  {
    cout.width(20);
    cout << std::left << user_ctx.input.regions.names[i] << ": { ";
    for (unsigned j = 0; j < user_ctx.input.regions.tags[i].size(); ++j)
    {
      cout << user_ctx.input.regions.tags[i][j];
      if (j != user_ctx.input.regions.tags[i].size()-1)
        cout  << ", ";
      else
        cout << " }" << endl;
    }
  }
  cout << endl;

  cout << "Parameters:\n";
  cout << "-----------\n";
  for (unsigned i = 0; i < user_ctx.input.parameters.size(); ++i)
  {
    cout.width(20);
    cout << std::left << user_ctx.input.parameters[i].name << ": { ";
    for (unsigned j = 0; j < user_ctx.input.parameters[i].values.size() ; ++j)
    {
      cout << user_ctx.input.parameters[i].values[j];
      if (j != user_ctx.input.parameters[i].values.size()-1)
        cout << ",  ";
      else
        cout << " }" << endl;
    }
  }


}

int main(int argc, char **argv)
{

  AppCtx user_ctx;

  if (argc > 1)
    user_ctx.input.cfg_fname = argv[1];
  else
  {
    printf("\nUsage: %s <path_to_your_cfg_file.lua> \n\nExample:\n%s some_dir/test.lua\n\n", argv[0], argv[0]);
    return 1;
  }

  if( read_input(user_ctx.input) )
    throw;

  int    petsc_argc = (int) user_ctx.input.petsc_opts.tokens.size();
  char** petsc_argv =       user_ctx.input.petsc_opts.tokens.data();

  PetscInitialize(&petsc_argc, &petsc_argv, PETSC_NULL, help);

  print_settings(user_ctx);
  
  user_ctx.mp = new MeshT(user_ctx.input.settings.space_dim);
  user_ctx.initAll();

  PetscFinalize();
  return 0;
}





























