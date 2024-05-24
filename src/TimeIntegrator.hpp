#ifndef TIME_INTEGRATOR_HPP_
#define TIME_INTEGRATOR_HPP_

#include <vector>

#include "Enums.hpp"
#include "ModularityUtils.hpp"
#include "ConservativeVariablesCC.hpp"
#include "GodunovFlux.hpp"
#include "ConstainedTransport.hpp"
#include "AddGhostCells.hpp"

ConservativeVariablesCC EulerAdvance(const ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int order, int nghost, Riemann rs);

ConservativeVariablesCC Euler(ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int order, int nghost, Riemann rs);

ConservativeVariablesCC TVDRK2(ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int order, int nghost, Riemann rs);

ConservativeVariablesCC TVDRK3(ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int order, int nghost, Riemann rs);

#endif // TIME_INTEGRATOR_HPP_