#ifndef CONSERVATIVE_VARIABLES_CC_HPP_
#define CONSERVATIVE_VARIABLES_CC_HPP_

#include <vector>
#include <string>

#include "PrimitiveVariablesCC.hpp"

class ConservativeVariablesCC {
public:
    int nx;
    int ny;
    std::vector<std::vector<double>> rho;
    std::vector<std::vector<double>> rhovx;
    std::vector<std::vector<double>> rhovy;
    std::vector<std::vector<double>> rhovz;
    std::vector<std::vector<double>> Bx;
    std::vector<std::vector<double>> By;
    std::vector<std::vector<double>> Bz;
    const double P = 10.0;

    // Constructor that takes grid dimensions as arguments
    ConservativeVariablesCC(const PrimitiveVariablesCC& P_cc) : nx(P_cc.nx), ny(P_cc.ny)
    {
        rho.resize(nx, std::vector<double>(ny));
        rhovx.resize(nx, std::vector<double>(ny));
        rhovy.resize(nx, std::vector<double>(ny));
        rhovz.resize(nx, std::vector<double>(ny));
        Bx.resize(nx, std::vector<double>(ny));
        By.resize(nx, std::vector<double>(ny));
        Bz.resize(nx, std::vector<double>(ny));

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                rho[i][j] = P_cc.rho[i][j];
                rhovx[i][j] = P_cc.rho[i][j] * P_cc.vx[i][j];
                rhovy[i][j] = P_cc.rho[i][j] * P_cc.vy[i][j];
                rhovz[i][j] = P_cc.rho[i][j] * P_cc.vz[i][j];
                Bx[i][j] = P_cc.Bx[i][j];
                By[i][j] = P_cc.By[i][j];
                Bz[i][j] = P_cc.Bz[i][j];
            }
        }
    }

    ~ConservativeVariablesCC() = default;
};

#endif // CONSERVATIVE_VARIABLES_CC_HPP_
