#include "GodunovFlux.hpp"

std::vector<ReconstructedValues> GodunovFluxX(const PrimitiveVariables& P_cc, int nghost, Reconstruction rec, Slope sl, Riemann rs){
    RiemannSolverFunction ChosenRiemannSolver = getRiemannSolver(rs);
    std::vector<ReconstructedValues> NumFluxx;

    for(int j = nghost; j < P_cc.ny - nghost; j++){
        for(int i = nghost; i < P_cc.nx + 1 - nghost; i++){
            NumFluxx.push_back(ChosenRiemannSolver(Interface(P_cc, i, j, rec, sl, nghost, Dir::X)));
        }
    }
    return NumFluxx;
}

std::vector<ReconstructedValues> GodunovFluxY(const PrimitiveVariables& P_cc, int nghost, Reconstruction rec, Slope sl, Riemann rs){
    RiemannSolverFunction ChosenRiemannSolver = getRiemannSolver(rs);
    std::vector<ReconstructedValues> NumFluxy;
    
    for(int i = nghost; i < P_cc.nx - nghost; i++){
        for(int j = nghost; j < P_cc.ny + 1 - nghost; j++){
            NumFluxy.push_back(ChosenRiemannSolver(Interface(P_cc, i, j, rec, sl, nghost, Dir::Y)));
        }
    }
    return NumFluxy;
}

ConservativeVariables ComputeFluxDifferenceX(PrimitiveVariables& P_cc, int nghost, Reconstruction rec, Slope sl, Riemann rs){
    std::vector<ReconstructedValues> NumFluxx = GodunovFluxX(P_cc, nghost, rec, sl, rs);

    ConservativeVariables FluxDifx(P_cc.nx, P_cc.ny);

    for(int j = 0; j < P_cc.ny - 2*nghost; j++){
        for(int i = 0; i < P_cc.nx - 2*nghost; i++){
            FluxDifx.setflux(NumFluxx[(i+1)+(P_cc.nx + 1 - 2*nghost)*j] - NumFluxx[i+(P_cc.nx + 1 - 2*nghost)*j], i + nghost, j + nghost);
        }
    }
    return FluxDifx;
}

ConservativeVariables ComputeFluxDifferenceY(PrimitiveVariables& P_cc, int nghost, Reconstruction rec, Slope sl, Riemann rs){
    std::vector<ReconstructedValues> NumFluxy = GodunovFluxY(P_cc, nghost, rec, sl, rs);

    ConservativeVariables FluxDify(P_cc.nx, P_cc.ny);

    for(int i = 0; i < P_cc.nx - 2*nghost; i++){
        for(int j=0; j < P_cc.ny - 2*nghost; j++){
            FluxDify.setflux(NumFluxy[(j+1)+(P_cc.ny + 1 - 2*nghost)*i] - NumFluxy[j+(P_cc.ny + 1 - 2*nghost)*i], i + nghost, j + nghost);
        }
    }

    return FluxDify;
}