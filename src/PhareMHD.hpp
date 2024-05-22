#ifndef PHARE_MHD_HPP_
#define PHARE_MHD_HPP_

#include <iostream>
#include <string>
#include <sstream>

#include "PrimitiveVariablesCC.hpp"
#include "AddGhostCells.hpp"
#include "ComputeNewDt.hpp"
#include "TimeIntegrator.hpp"
#include "WrittingUtils.hpp"

enum Reconstruction { Constant };
enum RiemannSolver { Rusanov };
enum CTMethod { Average };
enum Integrator { EulerIntegrator, TVDRK2Integrator, TVDRK3Integrator };

// Function pointer typedef for the integrators
typedef ConservativeVariablesCC (*IntegratorFunction)(ConservativeVariablesCC&, double, double, double, int, int);
IntegratorFunction getIntegrator(Integrator intg);

void PhareMHD(const PrimitiveVariablesCC& P0cc, std::string resultDir, int order, int nghost, Reconstruction rec, RiemannSolver rs, CTMethod ct, Integrator intg, double Dx, double Dy, double FinalTime, double Dt = 0.0);

#endif // PHARE_MHD_HPP_
