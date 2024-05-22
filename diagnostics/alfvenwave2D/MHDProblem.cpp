#include "Initialisation.hpp"
#include "PhareMHD.hpp"

/*#include "PrimitiveVariablesCC.hpp"
#include "AddGhostCells.hpp"
#include "ComputeNewDt.hpp"
#include "TimeIntegrator.hpp"
#include "WrittingUtils.hpp"*/


#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

int main(){
    Initialisation I;

    std::string resultsDir = "results/";

    system(("rm -rf " + resultsDir).c_str());
    system(("mkdir -p " + resultsDir).c_str());

    PrimitiveVariablesCC P0cc(I.nx, I.ny);
    P0cc.init(I.rho, I.vx, I.vy, I.vz, I.Bx, I.By, I.Bz, I.P); //OK

    PhareMHD(P0cc, resultsDir, I.order, I.nghost, Constant, Rusanov, Average, EulerIntegrator, I.Dx, I.Dy, I.FinalTime, I.Dt);


    /*PrimitiveVariablesCC P0 = InitialiseGhostCells(P0cc, I.nghost);

    // Save initial values
    savePrimitiveVariables(P0, resultsDir + "Pcc_0.txt", I.nghost);

    ConservativeVariablesCC U0(P0); 
    saveConcervativeVariables(U0, resultsDir + "URK2_0_0.txt", I.nghost);
    UpdateGhostCells(U0, I.nghost);

    ConservativeVariablesCC Un1 = U0;

    double time = 0.0;

    if(I.Dt == 0.0){
        double Dt = ComPuteNewDt(U0, I.Dx, I.Dy, I.nghost);;
        int step = 1;

        while(time <= I.FinalTime){
            Un1 = TVDRK2(Un1, I.Dx, I.Dy, Dt, I.order, I.nghost);
            //UpdateGhostCells(Un1, I.nghost);

            time = time + Dt;
            Dt = ComPuteNewDt(Un1, I.Dx, I.Dy, I.nghost);
            std::cout<<time<<" "<<Dt<<std::endl;

            std::ostringstream filename;
            filename << resultsDir << "URK2_" << formatTime(time) << ".txt";
            saveConcervativeVariables(Un1, filename.str(), I.nghost);
        }
    }else{
        for(int step = 1; step * I.Dt <= I.FinalTime; step++){
            Un1 = TVDRK2(Un1, I.Dx, I.Dy, I.Dt, I.order, I.nghost);

            time = time + I.Dt;
            std::cout<<time<<std::endl;

            std::ostringstream filename;
            filename << resultsDir << "URK2_" << step << ".txt";
            saveConcervativeVariables(Un1, filename.str(), I.nghost);
        }
    }*/
}