#ifndef WRITTING_UTILS_HPP_
#define WRITTING_UTILS_HPP_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "ConservativeVariables.hpp"
#include "PrimitiveVariables.hpp"

inline void savePrimitiveVariables(const PrimitiveVariables& P_cc, const std::string& filename, int nghost) {
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

inline void saveConcervativeVariables(const ConservativeVariables& C_cc, const std::string& filename, int nghost) {
    std::ofstream outFile(filename);

    if (outFile.is_open()) {
        for (int j = nghost; j < C_cc.ny - nghost; ++j) {
            for (int i = nghost; i < C_cc.nx - nghost; ++i) {
                outFile << std::setw(15) << C_cc.rho[j][i] << " "
                        << std::setw(15) << C_cc.rhovx[j][i] << " "
                        << std::setw(15) << C_cc.rhovy[j][i] << " "
                        << std::setw(15) << C_cc.rhovz[j][i] << " "
                        << std::setw(15) << C_cc.Bx[j][i] << " "
                        << std::setw(15) << C_cc.By[j][i] << " "
                        << std::setw(15) << C_cc.Bz[j][i] << " "
                        << std::setw(15) << C_cc.Etot[j][i] << std::endl;
            }
        }
        outFile.close();
        std::cout << "Wrote file: " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

inline void saveVectorToFile(const std::vector<ReconstructedValues>& values, const std::string& filename) {
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

inline void saveVectorToFile(const std::vector<double>& values, const std::string& filename) {
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

inline void saveVectorToFile(const std::vector<std::vector<double>>& values, const std::string& filename) {
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

inline std::string formatTime(double time) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(10) << time; // Set precision to include the decimal part
    std::string str = oss.str();
    
    // Replace the dot with an underscore
    std::replace(str.begin(), str.end(), '.', '_');
    
    return str;
}

#endif //WRITTING_UTILS_HPP_