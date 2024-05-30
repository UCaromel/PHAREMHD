#ifndef MODULARITY_UTILS_HPP_
#define MODULARITY_UTILS_HPP_

#include "Enums.hpp"
#include "TimeIntegrator.hpp"
#include "RiemannSolver.hpp"
#include "ConstainedTransport.hpp"

// Function pointer typedef for the riemann solvers
typedef ReconstructedValues (*RiemannSolverFunction)(const Interface&);
RiemannSolverFunction getRiemannSolver(Riemann rs);

// Function pointer typedef for the CT
typedef std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> (*CTFunction)(const ConservativeVariablesCC&, double, double, double, int, Reconstruction);
CTFunction getCT(CTMethod ct);

// Function pointer typedef for the integrators
typedef ConservativeVariablesCC (*IntegratorFunction)(ConservativeVariablesCC&, double, double, double, int, Reconstruction, Riemann, CTMethod);
IntegratorFunction getIntegrator(Integrator intg);

#endif