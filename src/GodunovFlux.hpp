#ifndef GODUNOV_FLUX_HPP_
#define GODUNOV_FLUX_HPP_

#include <iostream>
#include <vector>

#include "ModularityUtils.hpp"
#include "ReconstructedValues.hpp"
#include "Interface.hpp"
#include "RiemannSolver.hpp"

std::vector<ReconstructedValues> GodunovFluxX(const PrimitiveVariablesCC& P_cc, int order, int nghost, Riemann rs);

std::vector<ReconstructedValues> GodunovFluxY(const PrimitiveVariablesCC& P_cc, int order, int nghost, Riemann rs);

ConservativeVariablesCC ComputeFluxDifferenceX(PrimitiveVariablesCC& P_cc, int order, int nghost, Riemann rs);

ConservativeVariablesCC ComputeFluxDifferenceY(PrimitiveVariablesCC& P_cc, int order, int nghost, Riemann rs);

#endif //GODUNOV_FLUX_HPP_