#include "PhareMHD.hpp"

void PhareMHD(const PrimitiveVariablesCC& P0cc, std::string resultDir, int order, int nghost, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct, Integrator intg, double Dx, double Dy, double FinalTime, double Dt, int dumpfrequency) {
    PrimitiveVariablesCC P0 = InitialiseGhostCells(P0cc, nghost);
    ConservativeVariablesCC U0(P0);
    UpdateGhostCells(U0, nghost);

    ConservativeVariablesCC Un1 = U0;

    std::ostringstream filename0;
    filename0 << resultDir << "URK2_0_0.txt";
    saveConcervativeVariables(Un1, filename0.str(), nghost);

    double time = 0.0;
    IntegratorFunction ChosenIntegrator = getIntegrator(intg);

    if(Dt == 0.0){
        double Dt = ComPuteNewDt(U0, Dx, Dy, nghost);
        int step = 1;

        while(time <= FinalTime){
            Un1 = ChosenIntegrator(Un1, Dx, Dy, Dt, nghost, rec, sl, rs, ct);

            time = time + Dt;
            Dt = ComPuteNewDt(Un1, Dx, Dy, nghost);
            std::cout<<time<<" "<<Dt<<std::endl;

            if(step%dumpfrequency==0 || time >= FinalTime){
                std::ostringstream filename;
                filename << resultDir << "URK2_" << formatTime(time) << ".txt";
                saveConcervativeVariables(Un1, filename.str(), nghost);
            }
        }
    }else{
        for(int step = 1; step * Dt <= FinalTime; step++){
            Un1 = ChosenIntegrator(Un1, Dx, Dy, Dt, nghost, rec, sl, rs, ct);

            time = time + Dt;
            std::cout<<time<<std::endl;

            if(step%dumpfrequency==0 || time >= FinalTime){
                std::ostringstream filename;
                filename << resultDir << "URK2_" << step << ".txt";
                saveConcervativeVariables(Un1, filename.str(), nghost);
            }
        }
    }
}
