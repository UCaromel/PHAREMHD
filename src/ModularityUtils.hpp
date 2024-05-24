#ifndef MODULARITY_UTILS_HPP_
#define MODULARITY_UTILS_HPP_

#include "Enums.hpp"
#include "TimeIntegrator.hpp"
#include "RiemannSolver.hpp"

// Function pointer typedef for the riemann solvers
typedef ReconstructedValues (*RiemannSolverFunction)(const Interface&);
RiemannSolverFunction getRiemannSolver(Riemann rs);

// Function pointer typedef for the integrators
typedef ConservativeVariablesCC (*IntegratorFunction)(ConservativeVariablesCC&, double, double, double, int, int, Riemann);
IntegratorFunction getIntegrator(Integrator intg);

#endif