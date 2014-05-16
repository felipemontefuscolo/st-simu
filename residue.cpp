#include "common.hpp"



PetscErrorCode AppCtx::formJacobian(SNES snes, Vec Vec_u, Mat* /*Mat_jac*/, Mat* /*prejac*/, MatStructure* /*flag*/)
{
  PetscBool          found = PETSC_FALSE;
  char               snes_type[PETSC_MAX_PATH_LEN];
  
  PetscOptionsGetString(PETSC_NULL,"-snes_type",snes_type,PETSC_MAX_PATH_LEN-1,&found);
  
  if (found)
    if (string(snes_type) == string("test"))
    {
      cout << "WARNING: TESTING JACOBIAN !!!!! \n";
      this->formFunction(snes, Vec_u, systems[current_system_to_solve].Vec_res);
    }
  
  PetscFunctionReturn(0);

}



// form function of the system 0 (velocity and pressure)
PetscErrorCode formFunction_0(SNES /*snes*/, Vec /*x*/, Vec /*f*/)
{
  
  PetscFunctionReturn(0);
}

// form function of the system 0 (elasticity)
PetscErrorCode formFunction_1(SNES /*snes*/, Vec /*x*/, Vec /*f*/)
{
  
  PetscFunctionReturn(0);
}


PetscErrorCode AppCtx::formFunction(SNES snes, Vec x, Vec f)
{
  switch (current_system_to_solve)
  {
    case 0: // solve system 0
      formFunction_0(snes, x, f);
      break;
    case 1:
      formFunction_1(snes, x, f);
      break;
    default:
      throw;
  }  
  
  PetscFunctionReturn(0);
}



