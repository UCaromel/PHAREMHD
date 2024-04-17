#ifndef GODUNOV_FLUX_HPP_
#define GODUNOV_FLUX_HPP_

#include <vector>

#include "ReconstructedValues.hpp"
#include "Interface.hpp"
#include "RusanovRiemannSolver.hpp"

std::vector<ReconstructedValues> GodunovFluxx(PrimitiveVariablesCC& P_cc, int order){
    std::vector<ReconstructedValues> NumFluxx;
    for(int j=0; j<P_cc.ny; j++){
        for(int i=0; i<P_cc.nx+1; i++){
            NumFluxx.push_back(RusanovRiemannSolver(Interface(P_cc, i, j, order, Dir::X)));
        }
    }
    return NumFluxx;
}

std::vector<ReconstructedValues> GodunovFluxy(PrimitiveVariablesCC& P_cc, int order){
    std::vector<ReconstructedValues> NumFluxy;
    for(int i=0; i<P_cc.nx; i++){
        for(int j=0; j<P_cc.ny+1; j++){
            NumFluxy.push_back(RusanovRiemannSolver(Interface(P_cc, i, j, order, Dir::Y)));
        }
    }
    return NumFluxy;
}

ConservativeVariablesCC ComputeFluxDifference(PrimitiveVariablesCC& P_cc, int order){
    std::vector<ReconstructedValues> NumFluxx = GodunovFluxx(P_cc, order);
    std::vector<ReconstructedValues> NumFluxy = GodunovFluxy(P_cc, order);

    ConservativeVariablesCC FluxDifx(P_cc.nx, P_cc.ny);
    ConservativeVariablesCC FluxDify(P_cc.nx, P_cc.ny);
    ConservativeVariablesCC FluxDif(P_cc.nx, P_cc.ny);

    for(int j=0; j<P_cc.ny; j++){
        for(int i=0; i<P_cc.nx; i++){
            FluxDifx.set(NumFluxx[(i+1)+P_cc.nx*j] - NumFluxx[i+P_cc.nx*j], i, j);
        }
    }

    for(int i=0; i<P_cc.nx; i++){
        for(int j=0; j<P_cc.ny; j++){
            FluxDify.set(NumFluxy[(j+1)+P_cc.ny*i] - NumFluxy[j+P_cc.ny*i], i, j);
        }
    }

    FluxDif = FluxDifx - FluxDify;

    return FluxDif;
}


#endif //GODUNOV_FLUX_HPP_