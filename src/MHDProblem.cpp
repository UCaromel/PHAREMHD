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

void writeVectorToFile(const std::vector<ReconstructedValues>& values, const std::string& filename) {
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

    std::cout<<U0.rho[0][0]<<" "<<U0.rhovx[0][0]<<" "<<U0.rhovy[0][0]<<" "<<U0.rhovz[0][0]<<" "<<U0.Bx[0][0]<<" "<<U0.By[0][0]<<" "<<U0.Bz[0][0]<<std::endl; //OK
    std::cout<<U0(0,0).rho<<" "<<U0(0,0).vx<<" "<<U0(0,0).vy<<" "<<U0(0,0).vz<<" "<<U0(0,0).Bx<<" "<<U0(0,0).By<<" "<<U0(0,0).Bz<<std::endl; //OK

    ConservativeVariablesCC Uset(I.nx, I.ny); Uset.set(U0(0,0), 0, 0);
    std::cout<<Uset.rho[0][0]<<" "<<Uset.rhovx[0][0]<<" "<<Uset.rhovy[0][0]<<" "<<Uset.rhovz[0][0]<<" "<<Uset.Bx[0][0]<<" "<<Uset.By[0][0]<<" "<<Uset.Bz[0][0]<<std::endl; //OK

    ConservativeVariablesCC Ughost = AddGhostCells(U0, I.nghost); //OK
    saveConcervativeVariables(Ughost, resultsDir + "Ughost.txt");

    std::vector<ReconstructedValues> godunovfluxx = GodunovFluxX(P_cc, I.order, I.nghost); //OK
    writeVectorToFile(godunovfluxx, resultsDir + "godunovfluxx.txt");

    std::vector<ReconstructedValues> godunovfluxy = GodunovFluxY(P_cc, I.order, I.nghost); //OK
    writeVectorToFile(godunovfluxy, resultsDir + "godunovfluxy.txt");

    ConservativeVariablesCC fluxx = ComputeFluxDifferenceX(P_cc, I.order, I.nghost); //OK
    saveConcervativeVariables(fluxx, resultsDir + "fluxx.txt");

    ConservativeVariablesCC fluxy = ComputeFluxDifferenceY(P_cc, I.order, I.nghost); //OK
    saveConcervativeVariables(fluxy, resultsDir + "fluxy.txt");

    ConservativeVariablesCC Ueuler = EulerAdvance(U0, I.Dx, I.Dy, I.Dt, I.order, I.nghost); //OK
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