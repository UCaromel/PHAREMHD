#include "Initialisation.hpp"
#include "ConservativeVariablesCC.hpp"
#include "PrimitiveVariablesCC.hpp"
#include "TimeIntegrator.hpp"
#include "AddGhostCells.hpp"

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

    PrimitiveVariablesCC Pcc(I.nx, I.ny);
    Pcc.init(I.rho, I.vx, I.vy, I.vz, I.Bx, I.By, I.Bz, I.P); //OK

    PrimitiveVariablesCC P_cc = InitialiseGhostCells(Pcc, I.nghost);

    std::string resultsDir = "results/";

    // Create the results directory if it doesn't exist
    system(("mkdir -p " + resultsDir).c_str());

    // Save initial values
    savePrimitiveVariables(P_cc, resultsDir + "P_cc_initial.txt", I.nghost);

    ConservativeVariablesCC U0(P_cc); //OK
    saveConcervativeVariables(U0, resultsDir + "URK2_0.txt", I.nghost);
    UpdateGhostCells(U0, I.nghost);

    std::cout<<U0.rho[0][0]<<" "<<U0.rhovx[0][0]<<" "<<U0.rhovy[0][0]<<" "<<U0.rhovz[0][0]<<" "<<U0.Bx[0][0]<<" "<<U0.By[0][0]<<" "<<U0.Bz[0][0]<<std::endl; //OK
    std::cout<<U0(0,0).rho<<" "<<U0(0,0).vx<<" "<<U0(0,0).vy<<" "<<U0(0,0).vz<<" "<<U0(0,0).Bx<<" "<<U0(0,0).By<<" "<<U0(0,0).Bz<<std::endl; //OK

    ConservativeVariablesCC Uset(I.nx, I.ny); Uset.set(U0(0,0), 0, 0);
    std::cout<<Uset.rho[0][0]<<" "<<Uset.rhovx[0][0]<<" "<<Uset.rhovy[0][0]<<" "<<Uset.rhovz[0][0]<<" "<<Uset.Bx[0][0]<<" "<<Uset.By[0][0]<<" "<<Uset.Bz[0][0]<<std::endl; //OK

    std::vector<ReconstructedValues> godunovfluxx = GodunovFluxX(P_cc, I.order, I.nghost); //OK
    writeVectorToFile(godunovfluxx, resultsDir + "godunovfluxx.txt");

    std::vector<ReconstructedValues> godunovfluxy = GodunovFluxY(P_cc, I.order, I.nghost); //OK
    writeVectorToFile(godunovfluxy, resultsDir + "godunovfluxy.txt");

    ConservativeVariablesCC fluxx = ComputeFluxDifferenceX(P_cc, I.order, I.nghost); //OK
    saveConcervativeVariables(fluxx, resultsDir + "fluxx.txt", I.nghost);

    ConservativeVariablesCC fluxy = ComputeFluxDifferenceY(P_cc, I.order, I.nghost); //OK
    saveConcervativeVariables(fluxy, resultsDir + "fluxy.txt", I.nghost);

    ConservativeVariablesCC UnoCT= EulerAdvance(U0, I.Dx, I.Dy, I.Dt, I.order, I.nghost); //OK
    saveConcervativeVariables(UnoCT, resultsDir + "UnoCT.txt", I.nghost);

    ConservativeVariablesCC Ueuler = Euler(U0, I.Dx, I.Dy, I.Dt, I.order, I.nghost); //OK
    saveConcervativeVariables(Ueuler, resultsDir + "Ueuler.txt", I.nghost);

    ConservativeVariablesCC URK2 = TVDRK2(U0, I.Dx, I.Dy, I.Dt, I.order, I.nghost); //OK
    saveConcervativeVariables(URK2, resultsDir + "URK2.txt", I.nghost);

    ConservativeVariablesCC URK3 = TVDRK3(U0, I.Dx, I.Dy, I.Dt, I.order, I.nghost); //OK
    saveConcervativeVariables(URK3, resultsDir + "URK3.txt", I.nghost);

    // Likely a problem with the equation of state.
    PrimitiveVariablesCC PRK3(URK3);
    savePrimitiveVariables(PRK3, resultsDir + "PRK3.txt", I.nghost);

    std::cout<<"P : "<<PRK3.P[0][10]<<", "<<"P eos : "<<EosP(URK3(10,0))<<std::endl;

    for(int step = 1; step * I.Dt <= I.FinalTime; step++){
        ConservativeVariablesCC Un1 = TVDRK2(U0, I.Dx, I.Dy, I.Dt, I.order, I.nghost);

        std::ostringstream filename;
        filename << resultsDir << "URK2_" << step << ".txt";
        saveConcervativeVariables(Un1, filename.str(), I.nghost);
        U0 = Un1;
    }
}