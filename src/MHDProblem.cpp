#include "Initialisation.hpp"
#include "ConservativeVariablesCC.hpp"
#include "PrimitiveVariablesCC.hpp"
#include "TimeIntegrator.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>


void savePrimitiveVariables(const PrimitiveVariablesCC& P_cc, const std::string& filename) {
    std::ofstream outFile(filename);

    if (outFile.is_open()) {
        for (int i = 0; i < P_cc.nx; ++i) {
            for (int j = 0; j < P_cc.ny; ++j) {
                outFile << std::setw(15) << P_cc.rho[i][j] << " "
                        << std::setw(15) << P_cc.vx[i][j] << " "
                        << std::setw(15) << P_cc.vy[i][j] << " "
                        << std::setw(15) << P_cc.vz[i][j] << " "
                        << std::setw(15) << P_cc.Bx[i][j] << " "
                        << std::setw(15) << P_cc.By[i][j] << " "
                        << std::setw(15) << P_cc.Bz[i][j] << " "
                        << std::setw(15) << P_cc.P << std::endl;
            }
        }
        outFile.close();
        std::cout << "Wrote file: " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void saveConcervativeVariables(const ConservativeVariablesCC& P_cc, const std::string& filename) {
    std::ofstream outFile(filename);

    if (outFile.is_open()) {
        for (int i = 0; i < P_cc.nx; ++i) {
            for (int j = 0; j < P_cc.ny; ++j) {
                outFile << std::setw(15) << P_cc.rho[i][j] << " "
                        << std::setw(15) << P_cc.rhovx[i][j] << " "
                        << std::setw(15) << P_cc.rhovy[i][j] << " "
                        << std::setw(15) << P_cc.rhovz[i][j] << " "
                        << std::setw(15) << P_cc.Bx[i][j] << " "
                        << std::setw(15) << P_cc.By[i][j] << " "
                        << std::setw(15) << P_cc.Bz[i][j] << " "
                        << std::setw(15) << P_cc.P << std::endl;
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

    PrimitiveVariablesCC P_cc(I.nx, I.ny);
    P_cc.init(I.rho, I.vx, I.vy, I.vz, I.Bx, I.By, I.Bz, I.P); //OK

    std::string resultsDir = "results/";

    // Create the results directory if it doesn't exist
    system(("mkdir -p " + resultsDir).c_str());

    // Save initial values
    savePrimitiveVariables(P_cc, resultsDir + "P_cc_initial.txt");

    ConservativeVariablesCC U0(P_cc); //OK

    saveConcervativeVariables(U0, resultsDir + "U0.txt");

    std::cout<<U0.rho[0][0]<<" "<<U0.rhovx[0][0]<<" "<<U0.rhovy[0][0]<<" "<<U0.rhovz[0][0]<<" "<<U0.Bx[0][0]<<" "<<U0.By[0][0]<<" "<<U0.Bz[0][0]<<" "; //OK
    std::cout<<U0(0,0).rho<<" "<<U0(0,0).vx<<" "<<U0(0,0).vy<<" "<<U0(0,0).vz<<" "<<U0(0,0).Bx<<" "<<U0(0,0).By<<" "<<U0(0,0).Bz<<" "; //OK

    ConservativeVariablesCC Ughost = AddGhostCells(U0, I.nghost); //OK

    saveConcervativeVariables(Ughost, resultsDir + "Ughost.txt");

    

    // debugging to be done, uninitialised P somewhere.
    ConservativeVariablesCC Ueuler = EulerAdvance(U0, I.Dx, I.Dy, I.Dt, I.order, I.nghost);

    saveConcervativeVariables(Ueuler, resultsDir + "Ueuler.txt");

    /*for(int step = 0; step * I.Dt <= I.FinalTime; ++step){
        ConservativeVariablesCC U0(P_cc);

        // Write the vector components to files
        std::ostringstream filename;
        filename << resultsDir << "U0" << step << ".txt";
        savePrimitiveVariables(U0, filename.str());

        ConservativeVariablesCC Un1 = TVDRK2(U0, I.Dx, I.Dy, I.Dt, I.order, I.nghost);
        P_cc = PrimitiveVariablesCC(Un1);

        // Write the vector components to files
        std::ostringstream filename;
        filename << resultsDir << "P_cc_" << step << ".txt";
        savePrimitiveVariables(P_cc, filename.str());
    }*/
}