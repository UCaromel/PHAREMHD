#include "Initialisation.hpp"
#include "ConservativeVariablesCC.hpp"
#include "PrimitiveVariablesCC.hpp"
#include "TimeIntegrator.hpp"
#include "AddGhostCells.hpp"
#include "ComputeNewDt.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>


void savePrimitiveVariables(const PrimitiveVariablesCC& P_cc, const std::string& filename, int nghost) {
    std::ofstream outFile(filename);

    if (outFile.is_open()) {
        for (int j = nghost; j < P_cc.ny - nghost; ++j) {
            for (int i = nghost; i < P_cc.nx - nghost; ++i) {
                outFile << std::setw(15) << P_cc.rho[j][i] << " "
                        << std::setw(15) << P_cc.vx[j][i] << " "
                        << std::setw(15) << P_cc.vy[j][i] << " "
                        << std::setw(15) << P_cc.vz[j][i] << " "
                        << std::setw(15) << P_cc.Bx[j][i] << " "
                        << std::setw(15) << P_cc.By[j][i] << " "
                        << std::setw(15) << P_cc.Bz[j][i] << " "
                        << std::setw(15) << P_cc.P[j][i] << std::endl;
            }
        }
        outFile.close();
        std::cout << "Wrote file: " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void saveConcervativeVariables(const ConservativeVariablesCC& P_cc, const std::string& filename, int nghost) {
    std::ofstream outFile(filename);

    if (outFile.is_open()) {
        for (int j = nghost; j < P_cc.ny - nghost; ++j) {
            for (int i = nghost; i < P_cc.nx - nghost; ++i) {
                outFile << std::setw(15) << P_cc.rho[j][i] << " "
                        << std::setw(15) << P_cc.rhovx[j][i] << " "
                        << std::setw(15) << P_cc.rhovy[j][i] << " "
                        << std::setw(15) << P_cc.rhovz[j][i] << " "
                        << std::setw(15) << P_cc.Bx[j][i] << " "
                        << std::setw(15) << P_cc.By[j][i] << " "
                        << std::setw(15) << P_cc.Bz[j][i] << " "
                        << std::setw(15) << P_cc.Etot[j][i] << std::endl;
            }
        }
        outFile.close();
        std::cout << "Wrote file: " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

int main(){
    Initialisation I;

    std::string resultsDir = "space_results/";

    // Create the results directory if it doesn't exist
    system(("mkdir -p " + resultsDir).c_str());

    double Dt = 5e-4;
    double finalTime = 1;
    int stepDy = 0;

    int dumpFrequency = 10;
    double Dy = 0.02; int ny = 50;

    // Loop in x
    for(;(Dy > I.Dy / 32.0) && (ny < 1600); Dy /= 2.0, ny *= 2){
        double time = 0.0;

        int nx = 3;
        Initialisation Isim(nx, ny, I.Dx, Dy, Dt);
        std::cout<<"Dx : "<<Dy<<" nx : "<<ny<<std::endl;

        PrimitiveVariablesCC P0cc(nx, ny);
        P0cc.init(Isim.rho, Isim.vx, Isim.vy, Isim.vz, Isim.Bx, Isim.By, Isim.Bz, Isim.P); 
        PrimitiveVariablesCC P0 = InitialiseGhostCells(P0cc, I.nghost);
        ConservativeVariablesCC U0(P0);
        UpdateGhostCells(U0, I.nghost);
        ConservativeVariablesCC Un1 = U0;

        std::ostringstream f1;
        f1 << resultsDir << stepDy << "_" << "URK2_0.txt";
        saveConcervativeVariables(Un1, f1.str(), I.nghost);

        for(int step = 1; step * Dt <= finalTime; step++){
            Un1 = TVDRK2(Un1, I.Dx, Dy, Dt, I.order, I.nghost);

            time = time + Dt;
            std::cout<<time<<std::endl;

            if(step%dumpFrequency==0 || time >= finalTime){
                std::ostringstream filename;
                filename << resultsDir <<stepDy<<"_"<< "URK2_" << step << ".txt";
                saveConcervativeVariables(Un1, filename.str(), I.nghost);
            }
        }
        stepDy++;
    }
}