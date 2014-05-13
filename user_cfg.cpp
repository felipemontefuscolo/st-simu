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

// DEFINE HERE THE DOF MAPPERS
void defineSystems(AppCtx& ctx)
{
//  // do not change here ------------
//  UserInput const& input = ctx.input;
//  std::vector<UnkField> const& unk_fields = input.unk_fields;
//  std::vector<System>& systems = ctx.systems;
//  // -------------------------------
//
//
//  // user settings
//
//  systems.resize(2); // number of systems
//                     // one for u-p and other for mesh velocity and geometry
//
//  // attach mesh to the system ... usually is not necessary to change it.
//  for (unsigned i = 0; i < systems.size(); ++i)
//    systems[i].dof_mapper.attachMesh(ctx.mp);
//
//  // number of copies stored per time step
//  systems[0].n_copies_per_ts = 2;
//  systems[1].n_copies_per_ts = 4; // two for geometry and two for mesh_velocity
//
//  // defining the ordering of the degrees of freedom and creating systems
//
//  {
//    UnkField const& var = unk_fields[0]; // velocity
//    systems[0].dof_mapper.addVariable(var.name.c_str(), var.n_comps*ctx.shapes_c[0].numDofsPerVertex(),
//                                                        var.n_comps*ctx.shapes_c[0].numDofsInRidge(),
//                                                        var.n_comps*ctx.shapes_c[0].numDofsInFacet(),
//                                                        var.n_comps*ctx.shapes_c[0].numDofsInCell(),
//                                                        var.n_regions, var.n_tags.data(), var.tags.data());
//  }
//
//  {
//    UnkField const& var = unk_fields[1]; // pressure
//    systems[0].dof_mapper.addVariable(var.name.c_str(), var.n_comps*ctx.shapes_c[1].numDofsPerVertex(),
//                                                        var.n_comps*ctx.shapes_c[1].numDofsInRidge(),
//                                                        var.n_comps*ctx.shapes_c[1].numDofsInFacet(),
//                                                        var.n_comps*ctx.shapes_c[1].numDofsInCell(),
//                                                        var.n_regions, var.n_tags.data(), var.tags.data());
//  }
//
//
//  {
//    UnkField const& var = unk_fields[2]; // geometry and mesh velocity
//    systems[1].dof_mapper.addVariable(var.name.c_str(), var.n_comps*ctx.shapes_c[2].numDofsPerVertex(),
//                                                        var.n_comps*ctx.shapes_c[2].numDofsInRidge(),
//                                                        var.n_comps*ctx.shapes_c[2].numDofsInFacet(),
//                                                        var.n_comps*ctx.shapes_c[2].numDofsInCell(),
//                                                        var.n_regions, var.n_tags.data(), var.tags.data());
//  }
//
}
















