#ifndef CONSERVATIVE_VARIABLES_CC_HPP_
#define CONSERVATIVE_VARIABLES_CC_HPP_

#include <vector>
#include <string>

#include "PrimitiveVariablesCC.hpp"

class ConservativeVariablesCC {
public:
    int nx;
    int ny;
    std::vector<double> rho;
    std::vector<double> rhovx;
    std::vector<double> rhovy;
    std::vector<double> rhovz;
    std::vector<double> Bx;
    std::vector<double> By;
    std::vector<double> Bz;
    const double P = 10.0;

    // Constructor that takes grid dimensions as arguments
    ConservativeVariablesCC(const PrimitiveVariablesCC& P_cc) : nx(P_cc.nx), ny(P_cc.ny)
    {
        rho.resize(nx * ny);
        rhovx.resize(nx * ny);
        rhovy.resize(nx * ny);
        rhovz.resize(nx * ny);
        Bx.resize(nx * ny);
        By.resize(nx * ny);
        Bz.resize(nx * ny);

        for (int i = 0; i < nx * ny; i++) {
            rho[i] = P_cc.rho[i];
            rhovx[i] = P_cc.rho[i] * P_cc.vx[i];
            rhovy[i] = P_cc.rho[i] * P_cc.vy[i];
            rhovz[i] = P_cc.rho[i] * P_cc.vz[i];
            Bx[i] = P_cc.Bx[i];
            By[i] = P_cc.By[i];
            Bz[i] = P_cc.Bz[i];
        }
    }

    ~ConservativeVariablesCC() = default;
};



#endif // CONSERVATIVE_VARIABLES_CC_HPP_