#include "PhareMHD.hpp"

IntegratorFunction getIntegrator(Integrator intg) {
    switch (intg) {
        case EulerIntegrator: return &Euler;
        case TVDRK2Integrator: return &TVDRK2;
        case TVDRK3Integrator: return &TVDRK3;
        default: throw std::invalid_argument("Unknown integrator");
    }
}

void PhareMHD(const PrimitiveVariablesCC& P0cc, std::string resultDir, int order, int nghost, Reconstruction rec, RiemannSolver rs, CTMethod ct, Integrator intg, double Dx, double Dy, double FinalTime, double Dt) {
    PrimitiveVariablesCC P0 = InitialiseGhostCells(P0cc, nghost);
    ConservativeVariablesCC U0(P0);
    UpdateGhostCells(U0, nghost);

    ConservativeVariablesCC Un1 = U0;

    double time = 0.0;
    IntegratorFunction ChosenIntegrator = getIntegrator(intg);

    if(Dt == 0.0){
        double Dt = ComPuteNewDt(U0, Dx, Dy, nghost);;
        int step = 1;

        while(time <= FinalTime){
            Un1 = ChosenIntegrator(Un1, Dx, Dy, Dt, order, nghost);

            time = time + Dt;
            Dt = ComPuteNewDt(Un1, Dx, Dy, nghost);
            std::cout<<time<<" "<<Dt<<std::endl;

            std::ostringstream filename;
            filename << resultDir << "URK2_" << formatTime(time) << ".txt";
            saveConcervativeVariables(Un1, filename.str(), nghost);
        }
    }else{
        for(int step = 1; step * Dt <= FinalTime; step++){
            Un1 = ChosenIntegrator(Un1, Dx, Dy, Dt, order, nghost);

            time = time + Dt;
            std::cout<<time<<std::endl;

            std::ostringstream filename;
            filename << resultDir << "URK2_" << step << ".txt";
            saveConcervativeVariables(Un1, filename.str(), nghost);
        }
    }
}
