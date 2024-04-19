#ifndef GODUNOV_FLUX_HPP_
#define GODUNOV_FLUX_HPP_

#include <vector>

#include "ReconstructedValues.hpp"
#include "Interface.hpp"
#include "RusanovRiemannSolver.hpp"

std::vector<ReconstructedValues> GodunovFluxX(PrimitiveVariablesCC& P_cc, int order, int nghost){
    std::vector<ReconstructedValues> NumFluxx;
    PrimitiveVariablesCC Pghost = AddGhostCells(P_cc, nghost);
    for(int j=0; j<P_cc.ny; j++){
        for(int i=0; i<P_cc.nx+1; i++){
            NumFluxx.push_back(RusanovRiemannSolver(Interface(Pghost, i, j, order, nghost, Dir::X)));
        }
    }
    return NumFluxx;
}

std::vector<ReconstructedValues> GodunovFluxY(PrimitiveVariablesCC& P_cc, int order, int nghost){
    std::vector<ReconstructedValues> NumFluxy;
    PrimitiveVariablesCC Pghost = AddGhostCells(P_cc, nghost);
    for(int i=0; i<P_cc.nx; i++){
        for(int j=0; j<P_cc.ny+1; j++){
            NumFluxy.push_back(RusanovRiemannSolver(Interface(Pghost, i, j, order, nghost, Dir::Y)));
        }
    }
    return NumFluxy;
}

ConservativeVariablesCC ComputeFluxDifferenceX(PrimitiveVariablesCC& P_cc, int order, int nghost){
    std::vector<ReconstructedValues> NumFluxx = GodunovFluxX(P_cc, order, nghost);

    ConservativeVariablesCC FluxDifx(P_cc.nx, P_cc.ny);

    for(int j=0; j<P_cc.ny; j++){
        for(int i=0; i<P_cc.nx; i++){
            FluxDifx.set(NumFluxx[(i+1)+P_cc.nx*j] - NumFluxx[i+P_cc.nx*j], i, j);
        }
    }

    return FluxDifx;
}

ConservativeVariablesCC ComputeFluxDifferenceY(PrimitiveVariablesCC& P_cc, int order, int nghost){
    std::vector<ReconstructedValues> NumFluxy = GodunovFluxY(P_cc, order, nghost);

    ConservativeVariablesCC FluxDify(P_cc.nx, P_cc.ny);

    for(int i=0; i<P_cc.nx; i++){
        for(int j=0; j<P_cc.ny; j++){
            FluxDify.set(NumFluxy[(j+1)+P_cc.ny*i] - NumFluxy[j+P_cc.ny*i], i, j);
        }
    }

    return FluxDify;
}

#endif //GODUNOV_FLUX_HPP_