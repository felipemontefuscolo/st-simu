#include "common.hpp"

// Define here the exact solution for all variables.
// This function will also be used for initial conditions.
void exactSolution(int var, Vector const& X, double t, int tag, double const* pars, double* result)
{
  double Re = pars[1]/pars[2];
  int dim = X.size();

  switch (var)
  {
    case 0: // variable 0 = velocity
      {
        result[0] = (1.2 - X[1])*(1.2 + X[1])/(1.2*1.2);
        result[1] = 0.0;
      }
      break;
    case 1: // pressure
      {
        result[0] = 2.*(1.2 - X[0])/(1.2*1.2*Re);
      }
      break;
    case 2: // mesh_velocity
      {
        for (int i = 0; i < dim; ++i)
          result[i] = 0.0;
      }
      break;
    case 3: // geometry
      {
        for (int i = 0; i < dim; ++i)
          result[i] = X[i];
      }
      break;
    default:
      return;
  }
}















