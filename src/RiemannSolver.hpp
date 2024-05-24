#ifndef RIEMANN_SOLVER_HPP_
#define RIEMANN_SOLVER_HPP_

#include "Interface.hpp"

inline ReconstructedValues RusanovRiemannSolver(const Interface& inter){
    return 0.5 * (inter.fL + inter.fR) - 0.5 * inter.Splus * (inter.uR - inter.uL);
}

inline ReconstructedValues HLLRiemannSolver(const Interface& inter){
    if(inter.SL >= 0){
        return inter.fL;
    }
    else if((inter.SL <= 0) && (inter.SR >= 0)){
        return (inter.SR * inter.fL - inter.SL * inter.fR + inter.SL * inter.SR * (inter.uR - inter.uL)) / (inter.SR - inter.SL);
    }
    else{
        return inter.fR;
    }
}


#endif //RIEMANN_SOLVER_HPP_