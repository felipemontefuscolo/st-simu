

settings = {

  name = "Contact Angle Navier-Stokes",

  space_dim = 2,

  -- continuation parameter
  times = {from=2, step=1, to=4},

  stop_at_steady_state = false,

  steady_state_tolerance = 1.e-7, -- defined as |u^n+1 - u^n|/|u^n|

  input_meshfile = "some_dir/couette.msh",

  output_basename = "some_dir/couette.vtk",

  -- results print step
  print_step = 2,

  degree_of_pol_integrated_exactly_for_cells   = 3,
  degree_of_pol_integrated_exactly_for_facets  = 2,
  degree_of_pol_integrated_exactly_for_ridges  = 1,
  degree_of_pol_integrated_exactly_for_errors  = 6,  -- for error computation

  print_errors = false,

  -- Does not include the errors. The user must compute these quantities inside the program. Quantities will be printed at each time step.
  integ_quantities_names = {"volume", "kin_energy", "pot_energy", "surf_energy", "visc_power", "cl_dis", "solid_dis", "cl_max"},
  print_integ_quantities = {true,      false,        false,        true,          false,        false,    false,       true},

}

unknown_fields = {

  --[[
      for each unknown field, define:
      *  name : A string
      *  n_components : A positive integer
      *  interpolation : A string. Currently it can only be "Lagrange"
      *  degree : An integer. This is the degree of the interpolating polynomial (when is the case).
      *  system_number : An integer starting from 0. It's the linear system id where field belongs. If you define a field in 0 and another in 2, there must exists a third in 1.
      *  where_is_defined: A table of tables with tags where the field is defined. Example: {{1,2,3},{4,5}} has two regions. Empty {{}} mean everywhere
      *  dirichlet_tags: A table with the tags where Dirichlet b.c. is applied
  --]]

  -- field number 0 - velocity
  {
    name = "velocity",
    n_components = settings.space_dim,
    interpolation = "Lagrange",
    degree = 2,
    where_is_defined = {{}},
    dirichlet_tags = { 11,22,33 }
  },

  -- field number 1 - pressure
  {
    name = "pressure",
    n_components = 1,
    interpolation = "Lagrange",
    degree = 1,
    where_is_defined = {{1,2,3,4},{6,7,8}},
    dirichlet_tags = {  }
  },

  -- field number 2
  {
    name = "mesh_velocity",
    n_components = settings.space_dim,
    interpolation = "Lagrange",
    degree = 2,
    where_is_defined = {{}},
    dirichlet_tags = { 11,22,33 }
  },

  -- field number 3
  {
    name = "geometry",
    n_components = settings.space_dim,
    interpolation = "Lagrange",
    degree = 2,
    where_is_defined = {{}},
    dirichlet_tags = { 11,22,33 }
  }


}


systems = {

  --[[
  for each system define
  {
    active           =, -- true: allocate non-linear petsc solvers (SNES); false: allocated only petsc vectors
    first_dof_number  =, -- usually you don't have to change this. On could define a number greater 0 and manage the equations by yourself
    variables         =, -- variables of this system
    n_copiers_per_ts  =, -- number of copies of the vector per time step.
    mapper            =, -- "create_its_own": create a mapper wich maps mesh to dofs; sys_num: reuse dofs of the system sys_num
  },
  --]]
  
  -- system 0
  {
    active = true,
    first_dof_number = 0,
    fields = {0,1}, -- velocity and pressure
    n_copies_per_ts = 2,
    mapper = "create_its_own"
  },
  
  -- system 1
  {
    active = true,
    first_dof_number = 0,
    fields = {2}, -- mesh velocity
    n_copies_per_ts = 2,
    mapper = "create_its_own"
  },

  -- system 3
  {
    active = false,
    first_dof_number = 0,
    fields = {3}, -- geometry
    n_copies_per_ts = 2,
    mapper = 1
  }

}

region_marks = {

  {
    name = "solid tags",
    tags = {5,3,2,1}
  },

  {
    name = "gas tags",
    tags = {9,8,7,6}
  },


  do_after_reading = [[

    select("nodes", function(x,y,z,tag) return x end,  add)
    settag("selected", 2)

  ]],

}

parameters = {


 --[[
      The simulation will run for each combination of parameters.
      For example, if the parameters A and B are
        A = 1,
        B = {2,5}
      then the program will run for {A=1,B=2} and for {A=1,B=5}


      Accepted formats:

      <parameter name> =  X ,
      <parameter name> = {X0, X1, X2, ..., Xf},
      <parameter name> = {from = < X0 >, to = < Xf >, step = < dX >}

      this last format includes the Xf in the range.

 --]]

  {
    name = "theta",
    values = 0.5
  },

  {
    name   = "rho",
    values = 666
  },

  {
    name   = "mu",
    values = {from = 2.1, step = 2.2, to = 5*2.2},
  },

  {
    name   = "gamma",
    values = {1e-3, 1e-2, 1e-1, 1e+0},
  },

  {
    name   = "beta",
    values = {1e-3}
  }

}

petsc_options = {

  [[
    -pc_type lu lu                    # "Preconditioner" (one of) none jacobi pbjacobi bjacobi sor lu shell mg eisenstat ilu icc cholesky asm ksp composite redundant nn mat fieldsplit galerkin exotic openmp asa cp bfbt lsc redistribute tfs (PCSetType)
    sub_pc_type ilu
    -pc_factor_mat_solver_package mumps # MUMPS
    -mat_mumps_icntl_7 2
    pc_factor_levels 2
    sub_pc_factor_levels 1            # <0>  ativar quando o resolutor não estiver convergindo
    pc_composite_type multiplicative  # one of multiplicative, additive, special
    pc_composite_pcs ilu,ilu
    -ksp_type preonly gmres preonly   # (one of) preonly bcgs gmres cg cgne nash stcg gltr richardson chebychev tcqmr ibcgs bcgsl cgs tfqmr cr lsqr qcg bicg fgmres minres symmlq lgmres lcd broyden gcr (KSPSetType)
    ksp_initial_guess_nonzero 1       # não usar com precond asm+lu, alias nao usar nunca (talvez um erro da versão dev????)
    ksp_gmres_restart 300
    pc_factor_shift_type NONZERO
    pc_factor_shift_amount 1.e-10
    pc_factor_mat_ordering_type natural    # natural nd 1wd rcm qmd rowlength flow (PCFactorSetMatOrderingType)
    -pc_factor_reuse_ordering 1
    pc_factor_nonzeros_along_diagonal 1.e-10
    pc_factor_diagonal_fill

    pc_factor_fill 3.22746e-06
    pc_factor_in_place
    -ksp_rtol 1e-12 #<1e-8>

    fp_trap 1           # stop on floating-point exceptions
    on_error_abort 1
    log_trace stdout
    malloc_debug 1
    snes_fd 0
    -snes_ls basic      # line search: basic, cubic, quadratic
    -snes_type ksponly test   # Nonlinear solver method (one of) ls tr test picard ksponly vi ngmres sorqn
    snes_picard_alpha 1.2
    snes_test_display 1 # compare f.e. jacobian with f.d. jacobian
    snes_monitor_cancel 0 # cancela monitoramento
    -snes_monitor stdout
    snes_converged_reason 1
    snes_stol 1.e-15  # <1e-08>: Stop if step length less than
    snes_rtol 1.e-12  # <1e-08>: Stop if decrease in function norm less than
    ksp_monitor stdout
    snes_max_it 50
    mat_no_inode 1

    ######THREADS  #ativar os 3 primeiros somente
    vec_type seqpthread
    mat_type seqaijpthread
    thread_sync_type LOCKFREE
    vec_threads 3
    mat_threads 3
    use_thread_pool main
    nthreads 3

]]

}




