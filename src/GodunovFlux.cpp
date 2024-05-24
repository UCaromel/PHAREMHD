#include "GodunovFlux.hpp"

std::vector<ReconstructedValues> GodunovFluxX(const PrimitiveVariablesCC& P_cc, int order, int nghost, Riemann rs){
    RiemannSolverFunction ChosenRiemannSolver = getRiemannSolver(rs);
    std::vector<ReconstructedValues> NumFluxx;

    for(int j = nghost; j < P_cc.ny - nghost; j++){
        for(int i = nghost; i < P_cc.nx + 1 - nghost; i++){
            NumFluxx.push_back(ChosenRiemannSolver(Interface(P_cc, i, j, order, nghost, Dir::X)));
        }
    }
    return NumFluxx;
}

std::vector<ReconstructedValues> GodunovFluxY(const PrimitiveVariablesCC& P_cc, int order, int nghost, Riemann rs){
    RiemannSolverFunction ChosenRiemannSolver = getRiemannSolver(rs);
    std::vector<ReconstructedValues> NumFluxy;
    
    for(int i = nghost; i < P_cc.nx - nghost; i++){
        for(int j = nghost; j<P_cc.ny + 1 - nghost; j++){
            NumFluxy.push_back(ChosenRiemannSolver(Interface(P_cc, i, j, order, nghost, Dir::Y)));
        }
    }
    return NumFluxy;
}

ConservativeVariablesCC ComputeFluxDifferenceX(PrimitiveVariablesCC& P_cc, int order, int nghost, Riemann rs){
    std::vector<ReconstructedValues> NumFluxx = GodunovFluxX(P_cc, order, nghost, rs);

    ConservativeVariablesCC FluxDifx(P_cc.nx, P_cc.ny);

    for(int j = 0; j < P_cc.ny - 2*nghost; j++){
        for(int i = 0; i < P_cc.nx - 2*nghost; i++){
            FluxDifx.set(NumFluxx[(i+1)+(P_cc.nx + 1 - 2*nghost)*j] - NumFluxx[i+(P_cc.nx + 1 - 2*nghost)*j], i + nghost, j + nghost);
        }
    }
    return FluxDifx;
}

ConservativeVariablesCC ComputeFluxDifferenceY(PrimitiveVariablesCC& P_cc, int order, int nghost, Riemann rs){
    std::vector<ReconstructedValues> NumFluxy = GodunovFluxY(P_cc, order, nghost, rs);

    ConservativeVariablesCC FluxDify(P_cc.nx, P_cc.ny);

    for(int i = 0; i < P_cc.nx - 2*nghost; i++){
        for(int j=0; j < P_cc.ny - 2*nghost; j++){
            FluxDify.set(NumFluxy[(j+1)+(P_cc.ny + 1 - 2*nghost)*i] - NumFluxy[j+(P_cc.ny + 1 - 2*nghost)*i], i + nghost, j + nghost);
        }
    }

    return FluxDify;
}