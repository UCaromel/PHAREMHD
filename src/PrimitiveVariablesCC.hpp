#ifndef PRIMITIVE_VARIABLES_CC_HPP_
#define PRIMITIVE_VARIABLES_CC_HPP_

#include <vector>
#include <string>

#include "ReconstructedValues.hpp"
#include "ConservativeVariablesCC.hpp"

class PrimitiveVariablesCC {
public:
    int nx;
    int ny;
    std::vector<std::vector<double>> rho;
    std::vector<std::vector<double>> vx;
    std::vector<std::vector<double>> vy;
    std::vector<std::vector<double>> vz;
    std::vector<std::vector<double>> Bx;
    std::vector<std::vector<double>> By;
    std::vector<std::vector<double>> Bz;
    const double P = 10.0;

    // Constructor that takes grid dimensions as arguments
    PrimitiveVariablesCC(int nx_, int ny_) : nx(nx_), ny(ny_)
    {
        rho.resize(nx, std::vector<double>(ny));
        vx.resize(nx, std::vector<double>(ny));
        vy.resize(nx, std::vector<double>(ny));
        vz.resize(nx, std::vector<double>(ny));
        Bx.resize(nx, std::vector<double>(ny));
        By.resize(nx, std::vector<double>(ny));
        Bz.resize(nx, std::vector<double>(ny));
    }

    PrimitiveVariablesCC(const ConservativeVariablesCC& P_cc) : nx(P_cc.nx), ny(P_cc.ny)
    {
        rho.resize(nx, std::vector<double>(ny));
        vx.resize(nx, std::vector<double>(ny));
        vy.resize(nx, std::vector<double>(ny));
        vz.resize(nx, std::vector<double>(ny));
        Bx.resize(nx, std::vector<double>(ny));
        By.resize(nx, std::vector<double>(ny));
        Bz.resize(nx, std::vector<double>(ny));

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                rho[i][j] = P_cc.rho[i][j];
                vx[i][j] = P_cc.rhovx[i][j] / P_cc.rho[i][j];
                vy[i][j] = P_cc.rhovy[i][j] / P_cc.rho[i][j];
                vz[i][j] = P_cc.rhovz[i][j] / P_cc.rho[i][j];
                Bx[i][j] = P_cc.Bx[i][j];
                By[i][j] = P_cc.By[i][j];
                Bz[i][j] = P_cc.Bz[i][j];
            }
        }
    }


    ~PrimitiveVariablesCC() = default;

    void set(const ReconstructedValues& rv, int i, int j){
        rho[i][j] = rv.rho;
        vx[i][j] = rv.vx;
        vy[i][j] = rv.vy;
        vz[i][j] = rv.vz;
        Bx[i][j] = rv.Bx;
        By[i][j] = rv.By;
        Bz[i][j] = rv.Bz;
    }

    void init(const std::vector<std::vector<double>>& rho_, const std::vector<std::vector<double>>& vx_, const std::vector<std::vector<double>>& vy_, const std::vector<std::vector<double>>& vz_, const std::vector<std::vector<double>>& Bx_, const std::vector<std::vector<double>>& By_, const std::vector<std::vector<double>>& Bz_){
        rho = rho_;
        vx = vx_;
        vy = vy_;
        vz = vz_;
        Bx = Bx_;
        By = By_;
        Bz = Bz_;
    }

ReconstructedValues operator()(int i, int j) const {
    if (i < 0 || i >= nx || j < 0 || j >= ny)
        throw("Index out of range");

    return ReconstructedValues{rho[i][j], vx[i][j], vy[i][j], vz[i][j], Bx[i][j], By[i][j], Bz[i][j], P};
}

};

#endif // PRIMITIVE_VARIABLES_CC_HPP_
