#ifndef TIME_INTEGRATOR_HPP_
#define TIME_INTEGRATOR_HPP_

#include <vector>

#include <iostream>

#include "Enums.hpp"
#include "ModularityUtils.hpp"
#include "ConservativeVariables.hpp"
#include "GodunovFlux.hpp"
#include "ConstainedTransport.hpp"
#include "AddGhostCells.hpp"

ConservativeVariables EulerAdvance(const ConservativeVariables& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs);

ConservativeVariables Euler(ConservativeVariables& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct);

ConservativeVariables TVDRK2(ConservativeVariables& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct);

ConservativeVariables TVDRK3(ConservativeVariables& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct);

#endif // TIME_INTEGRATOR_HPP_