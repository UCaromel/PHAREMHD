// for higher order reconstructions later. For now it does midpoint but we need it to do PLM.

/*

#ifndef PRIMITIVE_VARIABLES_FC_HPP_
#define PRIMITIVE_VARIABLES_FC_HPP_

#include <vector>
#include <string>

#include "ConservativeVariablesCC.hpp"


class PrimitiveVariablesFC {
public:
    int nx;
    int ny;
    std::vector<double> rho;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;
    std::vector<double> Bx;
    std::vector<double> By;
    std::vector<double> Bz;
    const double P = 10.0;


    PrimitiveVariablesFC(const ConservativeVariablesCC& C_cc) : nx(C_cc.nx), ny(C_cc.ny)
    {
        // we only need interior values for the fluxes
        rho.resize(nx * ny);
        vx.resize((nx - 1) * ny);
        vy.resize(nx * (ny - 1));
        vz.resize(nx * ny);
        Bx.resize((nx - 1) * ny);
        By.resize(nx * (ny - 1));
        Bz.resize(nx * ny);

        rho = C_cc.rho;
        for(int i=0; i<nx-1; ++i){
            for(int j=0; j<ny; ++j){
                vx[i * ny + j] = 0.5*((C_cc.rhovx[i * ny + j] / C_cc.rho[i * ny + j]) + (C_cc.rhovx[(i + 1) * ny + j] / C_cc.rho[(i + 1) * ny + j]));
                Bx[i * ny + j] = 0.5*(C_cc.Bx[i * ny + j]  + C_cc.Bx[(i + 1) * ny + j] );
            }
        }

        for(int i=0; i<nx; ++i){
            for(int j=0; j<ny-1; ++j){
                vy[i * ny + j] = 0.5*((C_cc.rhovy[i * ny + j] / C_cc.rho[i * ny + j]) + (C_cc.rhovy[i * ny + j + 1] / C_cc.rho[i * ny + j + 1]));
                By[i * ny + j] = 0.5*(C_cc.By[i * ny + j]  + C_cc.By[i * ny + j + 1] );
            }
        }

        for(int i=0; i<nx*ny; ++i){
            vz[i] = C_cc.rhovz[i] / C_cc.rho[i];
        }
        Bz=C_cc.Bz;
    }

    ~PrimitiveVariablesFC() = default;
};

#endif // PRIMITIVE_VARIABLES_FC_HPP_

*/