#ifndef CONSERVATIVE_VARIABLES_CC_HPP_
#define CONSERVATIVE_VARIABLES_CC_HPP_

#include <vector>
#include <string>

#include "PrimitiveVariablesCC.hpp"

enum struct Dir {
    X,
    Y
};

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

        for (int i = 0; i < nx * ny; ++i) {
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
    
    // should problably be in a reconstruct to faces class as it is meant to be computed with ur and ul
    std::vector<double> ComputeFluxVector(int index, Dir dir){
        std::vector<double> flux(7, 0.0);
        if(dir == Dir::X){
            flux[0] = rhovx[index];
            flux[1] = rhovx[index] * (rhovx[index]/rho[index]) + P - Bx[index] * Bx[index];
            flux[2] = rhovx[index] * (rhovy[index]/rho[index]) - Bx[index] * By[index];
            flux[3] = rhovx[index] * (rhovz[index]/rho[index]) - Bx[index] * Bz[index];
            flux[4] = 0.0;
            flux[5] = By[index] * (rhovx[index] / rho[index]) - (rhovy[index] / rho[index]) * Bx[index];
            flux[6] = Bz[index] * (rhovx[index] / rho[index]) - (rhovz[index] / rho[index]) * Bx[index];
        }
        else if(dir == Dir::Y){
            flux[0] = rhovy[index];
            flux[1] = rhovy[index] * (rhovx[index]/rho[index]) - By[index] * Bx[index];
            flux[2] = rhovy[index] * (rhovy[index]/rho[index]) + P - By[index] * By[index];
            flux[3] = rhovy[index] * (rhovz[index]/rho[index]) - By[index] * Bz[index];
            flux[4] = Bx[index] * (rhovy[index] / rho[index]) - (rhovx[index] / rho[index]) * By[index];
            flux[5] = 0.0;
            flux[6] = Bz[index] * (rhovy[index] / rho[index]) - (rhovz[index] / rho[index]) * By[index];
        }
        return flux;
    }
};



#endif // CONSERVATIVE_VARIABLES_CC_HPP_