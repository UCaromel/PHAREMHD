#ifndef RIEMANN_SOLVER_HPP_
#define RIEMANN_SOLVER_HPP_

#include "ReconstructedValues.hpp"
#include "Interface.hpp"

class Interface;

ReconstructedValues RusanovRiemannSolver(const Interface& inter);

ReconstructedValues HLLRiemannSolver(const Interface& inter);

#endif //RIEMANN_SOLVER_HPP_