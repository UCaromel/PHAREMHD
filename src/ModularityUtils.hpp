#ifndef MODULARITY_UTILS_HPP_
#define MODULARITY_UTILS_HPP_

#include "Enums.hpp"
#include "SlopeLimiter.hpp"
#include "RiemannSolver.hpp"
#include "ConstainedTransport.hpp"
#include "TimeIntegrator.hpp"

class Interface;

// Function pointer typedef for the slope limiter
typedef ReconstructedValues (*SLFunction)(const ReconstructedValues&, const ReconstructedValues&);
SLFunction getSlopeLimiter(Slope sl);

// Function pointer typedef for the riemann solvers
typedef ReconstructedValues (*RiemannSolverFunction)(const Interface&);
RiemannSolverFunction getRiemannSolver(Riemann rs);

// Function pointer typedef for the CT
typedef std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> (*CTFunction)(const ConservativeVariablesCC&, double, double, double, int, Reconstruction, Slope, Riemann);
CTFunction getCT(CTMethod ct);

// Function pointer typedef for the integrators
typedef ConservativeVariablesCC (*IntegratorFunction)(ConservativeVariablesCC&, double, double, double, int, BoundaryConditions, Reconstruction, Slope, Riemann, CTMethod);
IntegratorFunction getIntegrator(Integrator intg);

#endif