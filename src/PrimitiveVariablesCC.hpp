#ifndef PRIMITIVE_VARIABLES_CC_HPP_
#define PRIMITIVE_VARIABLES_CC_HPP_

#include <vector>
#include <string>

#include "ReconstructedValues.hpp"
#include "ConservativeVariablesCC.hpp"

class ConservativeVariablesCC;

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
    double P;

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

    PrimitiveVariablesCC(const ConservativeVariablesCC& C_cc) : nx(C_cc.nx), ny(C_cc.ny)
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
                rho[i][j] = C_cc.rho[i][j];
                vx[i][j] = C_cc.rhovx[i][j] / C_cc.rho[i][j];
                vy[i][j] = C_cc.rhovy[i][j] / C_cc.rho[i][j];
                vz[i][j] = C_cc.rhovz[i][j] / C_cc.rho[i][j];
                Bx[i][j] = C_cc.Bx[i][j];
                By[i][j] = C_cc.By[i][j];
                Bz[i][j] = C_cc.Bz[i][j];
            }
        }
        P = C_cc.P;
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

    void init(const std::vector<std::vector<double>>& rho_, const std::vector<std::vector<double>>& vx_, const std::vector<std::vector<double>>& vy_, const std::vector<std::vector<double>>& vz_, const std::vector<std::vector<double>>& Bx_, const std::vector<std::vector<double>>& By_, const std::vector<std::vector<double>>& Bz_, const double P_){
        rho = rho_;
        vx = vx_;
        vy = vy_;
        vz = vz_;
        Bx = Bx_;
        By = By_;
        Bz = Bz_;
        P = P_;
    }

ReconstructedValues operator()(int i, int j) const {
    if (i < 0 || i >= nx || j < 0 || j >= ny)
        throw("Index out of range");

    return ReconstructedValues{rho[i][j], vx[i][j], vy[i][j], vz[i][j], Bx[i][j], By[i][j], Bz[i][j], P};
}

};

#endif // PRIMITIVE_VARIABLES_CC_HPP_
