#include "PhareMHD.hpp"

void PhareMHD(const PrimitiveVariables& P0cc, std::string resultDir, int order, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct, Integrator intg, double Dx, double Dy, double FinalTime, double Dt, dumpVariables dv, int dumpfrequency) {
    double DivB;

    PrimitiveVariables P0 = InitialiseGhostCells(P0cc, nghost, bc);
    ConservativeVariables U0(P0);
    UpdateGhostCells(U0, nghost, bc);

    DivB = CheckDivB(U0, Dx, Dy, nghost);
    std::cout<<"DivB = "<<DivB<<std::endl;

    ConservativeVariables Un1 = U0;

    if((dv == dumpVariables::Conservative) || (dv == dumpVariables::Both)){
        std::ostringstream filename0;
        filename0 << resultDir << "URK2_0_0.txt";
        saveConcervativeVariables(Un1, filename0.str(), nghost);
    } else if((dv == dumpVariables::Primitive) || (dv == dumpVariables::Both)){
        PrimitiveVariables Pdump0(Un1);
        std::ostringstream filename0;
        filename0 << resultDir << "PRK2_0_0.txt";
        savePrimitiveVariables(Pdump0, filename0.str(), nghost);
    }  else {
        throw std::invalid_argument("No valid variables selected for the dump");
    }

    double time = 0.0;
    IntegratorFunction ChosenIntegrator = getIntegrator(intg);

    if(Dt == 0.0){
        double Dt = ComPuteNewDt(U0, Dx, Dy, nghost);
        int step = 1;

        while(time <= FinalTime){
            Un1 = ChosenIntegrator(Un1, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct);

            time = time + Dt;
            Dt = ComPuteNewDt(Un1, Dx, Dy, nghost);
            std::cout<<time<<" "<<Dt<<std::endl;

            UpdateGhostCells(Un1, nghost, bc);
            DivB = CheckDivB(Un1, Dx, Dy, nghost);
            std::cout<<"DivB = "<<DivB<<std::endl;

            if((dv == dumpVariables::Conservative) || (dv == dumpVariables::Both)){
                if(step%dumpfrequency==0 || time >= FinalTime){
                    std::ostringstream filename;
                    filename << resultDir << "URK2_" << formatTime(time) << ".txt";
                    saveConcervativeVariables(Un1, filename.str(), nghost);
                }
            } else if((dv == dumpVariables::Primitive) || (dv == dumpVariables::Both)){
                if(step%dumpfrequency==0 || time >= FinalTime){
                    PrimitiveVariables Pdump(Un1);
                    std::ostringstream filename;
                    filename << resultDir << "PRK2_" << formatTime(time) << ".txt";
                    savePrimitiveVariables(Pdump, filename.str(), nghost);
                }
            } else {
                throw std::invalid_argument("No valid variables selected for the dump");
            }
            step++;
        }
    }else{
        for(int step = 1; step * Dt <= FinalTime; step++){
            Un1 = ChosenIntegrator(Un1, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct);

            time = time + Dt;
            std::cout<<time<<std::endl;

            UpdateGhostCells(Un1, nghost, bc);
            DivB = CheckDivB(Un1, Dx, Dy, nghost);
            std::cout<<"DivB = "<<DivB<<std::endl;

            if((dv == dumpVariables::Conservative) || (dv == dumpVariables::Both)){
                if(step%dumpfrequency==0 || time >= FinalTime){
                    std::ostringstream filename;
                    filename << resultDir << "URK2_" << formatTime(time) << ".txt";
                    saveConcervativeVariables(Un1, filename.str(), nghost);
                }
            } else if((dv == dumpVariables::Primitive) || (dv == dumpVariables::Both)){
                if(step%dumpfrequency==0 || time >= FinalTime){
                    PrimitiveVariables Pdump(Un1);
                    std::ostringstream filename;
                    filename << resultDir << "PRK2_" << formatTime(time) << ".txt";
                    savePrimitiveVariables(Pdump, filename.str(), nghost);
                }
            } else {
                throw std::invalid_argument("No valid variables selected for the dump");
            }
        }
    }
}
