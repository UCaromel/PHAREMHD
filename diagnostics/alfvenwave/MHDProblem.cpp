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

    std::string resultsDir = "results/";

    // Create the results directory if it doesn't exist
    system(("mkdir -p " + resultsDir).c_str());

    // Save initial values
    savePrimitiveVariables(P0, resultsDir + "Pcc_0.txt", I.nghost);

    ConservativeVariablesCC U0(P0); 
    saveConcervativeVariables(U0, resultsDir + "URK2_0.txt", I.nghost);
    UpdateGhostCells(U0, I.nghost);

    ConservativeVariablesCC Un1 = U0;
    ConstrainedTransport CT(U0, I.Dx, I.Dy, I.Dt, I.nghost);
    
    saveVectorToFile(CT.bx, resultsDir + "CTbx");
    saveVectorToFile(CT.by, resultsDir + "CTby");

    saveVectorToFile(CT.vx, resultsDir + "CTvx");
    saveVectorToFile(CT.vy, resultsDir + "CTvy");
    saveVectorToFile(CT.Bx, resultsDir + "CTBx");
    saveVectorToFile(CT.By, resultsDir + "CTBy");
    saveVectorToFile(CT.Ez, resultsDir + "CTEz");

    saveVectorToFile(CT.BX, resultsDir + "CTBX");
    saveVectorToFile(CT.BY, resultsDir + "CTBY");

    ConservativeVariablesCC U1 = EulerAdvance(U0, I.Dx, I.Dy, I.Dt, I.order, I.nghost);
    UpdateGhostCells(U1, I.nghost);

    ConstrainedTransport CT2(U1, I.Dx, I.Dy, I.Dt, I.nghost);

    saveVectorToFile(CT2.bx, resultsDir + "CT2bx");
    saveVectorToFile(CT2.by, resultsDir + "CT2by");

    saveVectorToFile(CT2.vx, resultsDir + "CT2vx");
    saveVectorToFile(CT2.vy, resultsDir + "CT2vy");
    saveVectorToFile(CT2.Bx, resultsDir + "CT2Bx");
    saveVectorToFile(CT2.By, resultsDir + "CT2By");
    saveVectorToFile(CT2.Ez, resultsDir + "CT2Ez");

    saveVectorToFile(CT2.BX, resultsDir + "CT2BX");
    saveVectorToFile(CT2.BY, resultsDir + "CT2BY");
/*
    double time = 0.0;

    if(I.Dt == 0.0){
        double Dt = ComPuteNewDt(U0, I.Dx, I.Dy, I.nghost);;
        int step = 1;

        while(time <= I.FinalTime){
            Un1 = Euler(U0, I.Dx, I.Dy, Dt, I.order, I.nghost);
            //UpdateGhostCells(Un1, I.nghost);

            time = time + Dt;
            Dt = ComPuteNewDt(Un1, I.Dx, I.Dy, I.nghost);
            std::cout<<time<<" "<<Dt<<std::endl;

            std::ostringstream filename;
            filename << resultsDir << "URK2_" << formatTime(time) << ".txt";
            saveConcervativeVariables(Un1, filename.str(), I.nghost);

            U0 = Un1;
        }
    }else{
        for(int step = 1; step * I.Dt <= I.FinalTime; step++){
            Un1 = TVDRK2(U0, I.Dx, I.Dy, I.Dt, I.order, I.nghost);
            time = time + I.Dt;
            std::cout<<time<<std::endl;

            std::ostringstream filename;
            filename << resultsDir << "URK2_" << step << ".txt";
            saveConcervativeVariables(Un1, filename.str(), I.nghost);

            U0 = Un1;
        }
    }*/
}