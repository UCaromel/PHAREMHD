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

void saveVectorToFile(const std::vector<ReconstructedValues>& values, const std::string& filename) {
    std::ofstream outFile(filename);

    if (outFile.is_open()) {
        for (const auto& value : values) {
            outFile << std::setw(15) << value.rho << " "
                    << std::setw(15) << value.vx << " "
                    << std::setw(15) << value.vy << " "
                    << std::setw(15) << value.vz << " "
                    << std::setw(15) << value.Bx << " "
                    << std::setw(15) << value.By << " "
                    << std::setw(15) << value.Bz << " "
                    << std::setw(15) << value.P << std::endl;
        }
        outFile.close();
        std::cout << "Wrote file: " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void saveVectorToFile(const std::vector<double>& values, const std::string& filename) {
    std::ofstream outFile(filename);

    if (outFile.is_open()) {
        for (const auto& value : values) {
            outFile << std::setw(15) << value <<std::endl;
        }
        outFile.close();
        std::cout << "Wrote file: " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void saveVectorToFile(const std::vector<std::vector<double>>& values, const std::string& filename) {
    std::ofstream outFile(filename);

    if (outFile.is_open()) {
        for (const auto& innerVector : values) {
            for (const auto& value : innerVector) {
                outFile << std::setw(15) << value;
            }
            outFile << std::endl;
        }
        outFile.close();
        std::cout << "Wrote file: " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

std::string formatTime(double time) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(10) << time; // Set precision to include the decimal part
    std::string str = oss.str();
    
    // Replace the dot with an underscore
    std::replace(str.begin(), str.end(), '.', '_');
    
    return str;
}

int main(){
    Initialisation I;

    PrimitiveVariablesCC P0cc(I.nx, I.ny);
    P0cc.init(I.rho, I.vx, I.vy, I.vz, I.Bx, I.By, I.Bz, I.P); //OK

    PrimitiveVariablesCC P0 = InitialiseGhostCells(P0cc, I.nghost);

    std::string resultsDir = "time_results/";

    // Create the results directory if it doesn't exist
    system(("mkdir -p " + resultsDir).c_str());

    ConservativeVariablesCC U0(P0); 
    UpdateGhostCells(U0, I.nghost);

    ConservativeVariablesCC UEuler = U0;
    ConservativeVariablesCC URK2 = U0;
    ConservativeVariablesCC URK3 = U0;

    int stepDt = 0;
    double finalTime = 0.001;

    int dumpFrequency = 10;

    for(double Dt = I.Dt; Dt >= I.Dt / 32.0; Dt /= 2.0){

        double time = 0.0;
        UEuler = U0;
        URK2 = U0;
        URK3 = U0;

        std::ostringstream f1;
        f1 << resultsDir << stepDt << "_" << "UEuler_0.txt";
        saveConcervativeVariables(UEuler, f1.str(), I.nghost);

        std::ostringstream f2;
        f2 << resultsDir << stepDt << "_" << "URK2_0.txt";
        saveConcervativeVariables(URK2, f2.str(), I.nghost);

        std::ostringstream f3;
        f3 << resultsDir << stepDt << "_" << "URK3_0.txt";
        saveConcervativeVariables(URK3, f3.str(), I.nghost);

        for(int step = 1; step * Dt <= finalTime; step++){
            UEuler = Euler(UEuler, I.Dx, I.Dy, Dt, I.order, I.nghost);
            URK2 = TVDRK2(URK2, I.Dx, I.Dy, Dt, I.order, I.nghost);
            URK3 = TVDRK3(URK3, I.Dx, I.Dy, Dt, I.order, I.nghost);

            time = time + Dt;
            std::cout<<time<<std::endl;

            if(step%dumpFrequency==0 || time >= finalTime){
                std::ostringstream filename;
                filename << resultsDir << stepDt << "_" << "UEuler_" << step << ".txt";
                saveConcervativeVariables(UEuler, filename.str(), I.nghost);

                std::ostringstream filename2;
                filename2 << resultsDir << stepDt << "_" << "URK2_" << step << ".txt";
                saveConcervativeVariables(URK2, filename2.str(), I.nghost);

                std::ostringstream filename3;
                filename3 << resultsDir << stepDt << "_" << "URK3_" << step << ".txt";
                saveConcervativeVariables(URK3, filename3.str(), I.nghost);
            }
        }
        stepDt++;
    }  
}