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

void AppCtx::destroyPetsc()
{
  // petsc objects
  //
  for (unsigned i = 0; i < systems.size(); ++i)
  {
    System& sys = systems[i];
    if (sys.active)
    {
      VecDestroy(&sys.Vec_res);
      MatDestroy(&sys.Mat_jac);
      SNESDestroy(&sys.snes);
    }
    for (unsigned j = 0; j < sys.Vec_u.size(); ++j)
      VecDestroy(&sys.Vec_u[j]);
  }
}

void AppCtx::printSettings() const
{
  //cout << "Cell type number " << CELL_TYPE << endl;
  //cout << "Cell type number " << FACET_TYPE << endl;
  //cout << "Cell type number " << RIDGE_TYPE << endl;

  cout << "\n\n*************************  " << settings.name << " *************************  " << endl << endl;

  cout << endl;
  cout << "Mesh data:" << endl;
  cout << "==========" << endl;
  cout << "input mesh file             : " << settings.input_meshfile << endl;
  cout << "output mesh file (basename) : " << settings.output_basename << endl;  
  cout << "# vertices   : " << mp->numVertices() << endl;
  cout << "# ridges(3D) : " << mp->numRidges() << endl;
  cout << "# facets     : " << mp->numFacets() << endl;
  cout << "# cells      : " << mp->numCells() << endl;

  cout << endl;
  cout << "Degrees of freedom:" << endl;
  cout << "===================" << endl;
  cout << "# dofs:" << endl;
  for (unsigned i = 0; i < systems.size(); ++i)
    if (systems[i].active)
      cout << "(active) system (" << i << ") : " << systems[i].n_dofs << endl;
    else
      cout << "         system (" << i << ") : " << systems[i].n_dofs << endl;

  cout << endl;
  cout << "Settings:\n";
  cout << "=========\n";
  cout << "mesh type                   : " << cell_name(CELL_TYPE) << endl;
  cout << "space dimension             : " << settings.space_dim << endl;
  cout << "time                        : " << "from=" << settings.time_from
                                           <<", to="<<settings.time_to
                                           <<", step=" << settings.time_step << endl;
  cout << "exact integ. degree (cell)  : " << settings.quadr_degree_cell << endl;
  cout << "exact integ. degree (facet) : " << settings.quadr_degree_facet << endl;
  cout << "exact integ. degree (ridge) : " << settings.quadr_degree_ridge << endl;

  cout << endl;

  cout << "Variables:\n";
  cout << "==========\n";
  {
    for (int j = 0; j < (int)unk_fields.size(); ++j)
    {
      AppCtx::UnkField const& var = unk_fields[j];
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
      if (var.where_is_defined.empty())
        cout << "{} }\n";

      cout << "dirichlet tags    : " << "{ ";
      for (unsigned i = 0; i < var.dirich_tags.size(); ++i)
      {
        cout << var.dirich_tags[i];
        if (i != var.dirich_tags.size()-1)
          cout << ",";
      }
      cout << " }\n" << endl;
    }
    cout << endl;
  }

  cout << "Regions:\n";
  cout << "========\n";
  for (unsigned i = 0; i < regions.names.size(); ++i)
  {
    cout.width(20);
    cout << std::left << regions.names[i] << ": { ";
    for (unsigned j = 0; j < regions.tags[i].size(); ++j)
    {
      cout << regions.tags[i][j];
      if (j != regions.tags[i].size()-1)
        cout  << ", ";
      else
        cout << " }" << endl;
    }
  }
  cout << endl;

  cout << "Parameters:\n";
  cout << "===========\n";
  for (unsigned i = 0; i < parameters.size(); ++i)
  {
    cout.width(20);
    cout << std::left << parameters[i].name << ": { ";
    for (unsigned j = 0; j < parameters[i].values.size() ; ++j)
    {
      cout << parameters[i].values[j];
      if (j != parameters[i].values.size()-1)
        cout << ",  ";
      else
        cout << " }" << endl;
    }
  }


}

void AppCtx::setUpInitialConditions()
{
  current_time = settings.time_from;
  
  // vertices dofs
  {
    VertexH v = mp->vertexBegin();
    VertexH v_end = mp->vertexEnd();
    for (; v != v_end; ++v)
    {
      if (v.isDisabled(mp))
        continue;
      
      
      
    }
    
  }
}

void AppCtx::solveTemporalProblem()
{
  setUpInitialConditions();
  
  
  
}

int main(int argc, char **argv)
{

  AppCtx user_ctx;

  if (argc > 1)
    user_ctx.cfg_fname = argv[1];
  else
  {
    printf("\nUsage: %s <path_to_your_cfg_file.lua> \n\nExample:\n%s some_dir/test.lua\n\n", argv[0], argv[0]);
    return 1;
  }

  if( read_input(user_ctx) )
    throw;

  int    petsc_argc = (int) user_ctx.petsc_opts.tokens.size();
  char** petsc_argv =       user_ctx.petsc_opts.tokens.data();

  PetscInitialize(&petsc_argc, &petsc_argv, PETSC_NULL, help);

  user_ctx.mp = new MeshT(user_ctx.settings.space_dim);
  user_ctx.initAll();
  user_ctx.printSettings();

  user_ctx.destroyPetsc();
  PetscFinalize();
  return 0;
}






























