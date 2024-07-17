#ifndef GODUNOV_FLUX_HPP_
#define GODUNOV_FLUX_HPP_

#include <iostream>
#include <vector>

#include "ModularityUtils.hpp"
#include "ReconstructedValues.hpp"
#include "Interface.hpp"
#include "RiemannSolver.hpp"

std::vector<ReconstructedValues> GodunovFluxX(const PrimitiveVariables& P_cc, double Dx, double Dy, int nghost, Reconstruction rec, Slope sl, Riemann rs, OptionalPhysics OptP);

std::vector<ReconstructedValues> GodunovFluxY(const PrimitiveVariables& P_cc, double Dx, double Dy, int nghost, Reconstruction rec, Slope sl, Riemann rs, OptionalPhysics OptP);

ConservativeVariables ComputeFluxDifferenceX(PrimitiveVariables& P_cc, double Dx, double Dy, int nghost, Reconstruction rec, Slope sl, Riemann rs, OptionalPhysics OptP);

ConservativeVariables ComputeFluxDifferenceY(PrimitiveVariables& P_cc, double Dx, double Dy, int nghost, Reconstruction rec, Slope sl, Riemann rs, OptionalPhysics OptP);

#endif //GODUNOV_FLUX_HPP_