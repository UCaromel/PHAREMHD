#ifndef PHARE_MHD_HPP_
#define PHARE_MHD_HPP_

#include <iostream>
#include <string>
#include <sstream>

#include "ModularityUtils.hpp"
#include "PrimitiveVariablesCC.hpp"
#include "AddGhostCells.hpp"
#include "ComputeNewDt.hpp"
#include "TimeIntegrator.hpp"
#include "WrittingUtils.hpp"


void PhareMHD(const PrimitiveVariablesCC& P0cc, std::string resultDir, int order, int nghost, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct, Integrator intg, double Dx, double Dy, double FinalTime, double Dt = 0.0, int dumpfrequency = 1);

#endif // PHARE_MHD_HPP_
