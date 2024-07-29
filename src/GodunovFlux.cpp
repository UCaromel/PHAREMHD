#include "GodunovFlux.hpp"

#include "WrittingUtils.hpp"

std::vector<ReconstructedValues> GodunovFluxX(const PrimitiveVariables& P_cc, double Dx, double Dy, int nghost, Reconstruction rec, Slope sl, Riemann rs, OptionalPhysics OptP){
    RiemannSolverFunction ChosenRiemannSolver = getRiemannSolver(rs);
    std::vector<ReconstructedValues> NumFluxx;

    for(int j = nghost; j < P_cc.ny - nghost; j++){
        for(int i = nghost; i < P_cc.nx + 1 - nghost; i++){
            auto inter = Interface(P_cc, i, j, Dx, Dy, nghost, rec, sl, OptP, Dir::X);
            auto riemann = ChosenRiemannSolver(inter);
            NumFluxx.push_back(riemann);
        }
    }
    return NumFluxx;
}

std::vector<ReconstructedValues> GodunovFluxY(const PrimitiveVariables& P_cc, double Dx, double Dy, int nghost, Reconstruction rec, Slope sl, Riemann rs, OptionalPhysics OptP){
    RiemannSolverFunction ChosenRiemannSolver = getRiemannSolver(rs);
    std::vector<ReconstructedValues> NumFluxy;
    
    for(int i = nghost; i < P_cc.nx - nghost; i++){
        for(int j = nghost; j < P_cc.ny + 1 - nghost; j++){
            NumFluxy.push_back(ChosenRiemannSolver(Interface(P_cc, i, j, Dx, Dy, nghost,  rec, sl, OptP, Dir::Y)));
        }
    }
    return NumFluxy;
}

ConservativeVariables ComputeFluxDifferenceX(PrimitiveVariables& P_cc, double Dx, double Dy, int nghost, Reconstruction rec, Slope sl, Riemann rs, OptionalPhysics OptP){
    std::vector<ReconstructedValues> NumFluxx = GodunovFluxX(P_cc, Dx, Dy, nghost, rec, sl, rs, OptP);

    ConservativeVariables FluxDifx(P_cc.nx, P_cc.ny);

    for(int j = 0; j < P_cc.ny - 2*nghost; j++){
        for(int i = 0; i < P_cc.nx - 2*nghost; i++){
            FluxDifx.set(NumFluxx[(i+1)+(P_cc.nx + 1 - 2*nghost)*j] - NumFluxx[i+(P_cc.nx + 1 - 2*nghost)*j], i + nghost, j + nghost);
        }
    }

    // static int i = 0;
    // std::ostringstream fdifx;
    // fdifx << "whislerwaveres/FluxDifx_" << i <<".txt";
    // saveConcervativeVariables(FluxDifx, fdifx.str(), nghost);
    // i++;

    return FluxDifx;
}

ConservativeVariables ComputeFluxDifferenceY(PrimitiveVariables& P_cc, double Dx, double Dy, int nghost, Reconstruction rec, Slope sl, Riemann rs, OptionalPhysics OptP){
    std::vector<ReconstructedValues> NumFluxy = GodunovFluxY(P_cc, Dx, Dy, nghost, rec, sl, rs, OptP);

    ConservativeVariables FluxDify(P_cc.nx, P_cc.ny);

    for(int i = 0; i < P_cc.nx - 2*nghost; i++){
        for(int j=0; j < P_cc.ny - 2*nghost; j++){
            FluxDify.set(NumFluxy[(j+1)+(P_cc.ny + 1 - 2*nghost)*i] - NumFluxy[j+(P_cc.ny + 1 - 2*nghost)*i], i + nghost, j + nghost);
        }
    }

    return FluxDify;
}