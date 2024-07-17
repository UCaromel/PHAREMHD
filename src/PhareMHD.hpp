#ifndef PHARE_MHD_HPP_
#define PHARE_MHD_HPP_

#include <iostream>
#include <string>
#include <sstream>

#include "Enums.hpp"
#include "ModularityUtils.hpp"
#include "PrimitiveVariables.hpp"
#include "AddGhostCells.hpp"
#include "ComputeNewDt.hpp"
#include "TimeIntegrator.hpp"
#include "ComputeJ.hpp"
#include "CheckDivB.hpp"
#include "WrittingUtils.hpp"
#include "PhysicalConstants.hpp"


void PhareMHD(const PrimitiveVariables& P0cc, std::string resultDir, int nghost,
              BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct, Integrator intg,
              double Dx, double Dy, double FinalTime, double Dt = 0.0, Consts cst = Consts(), OptionalPhysics OptP = Off,
              dumpVariables dv = Conservative, int dumpfrequency = 1);

#endif // PHARE_MHD_HPP_
