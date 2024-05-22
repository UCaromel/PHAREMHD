#ifndef RUSANOV_RIEMANN_SOLVER_HPP_
#define RUSANOV_RIEMANN_SOLVER_HPP_


#include "Interface.hpp"

inline ReconstructedValues RusanovRiemannSolver(const Interface& inter){
    return 0.5 * (inter.fL + inter.fR) - 0.5 * inter.Splus * (inter.uR - inter.uL);
}


#endif