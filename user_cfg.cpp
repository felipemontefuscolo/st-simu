#include "common.hpp"


void exactSolution(int var, Vector const& X, double t, double const* pars, double* result)
{

  switch (var)
  {
    case 0: // variable 0 = velocity
      {
        result[0] = pars[0]*X[0]*t;
      }
      break;
    default:
      throw;
  }

}














