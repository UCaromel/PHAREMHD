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
    std::vector<double> rho;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;
    std::vector<double> Bx;
    std::vector<double> By;
    std::vector<double> Bz;
    const double P = 10.0;

    // Constructor that takes grid dimensions as arguments
    PrimitiveVariablesCC(int nx_, int ny_) : nx(nx_), ny(ny_)
    {
        rho.resize(nx * ny);
        vx.resize(nx * ny);
        vy.resize(nx * ny);
        vz.resize(nx * ny);
        Bx.resize(nx * ny);
        By.resize(nx * ny);
        Bz.resize(nx * ny);
    }

    PrimitiveVariablesCC(const ConservativeVariablesCC& P_cc) : nx(P_cc.nx), ny(P_cc.ny)
    {
        rho.resize(nx * ny);
        vx.resize(nx * ny);
        vy.resize(nx * ny);
        vz.resize(nx * ny);
        Bx.resize(nx * ny);
        By.resize(nx * ny);
        Bz.resize(nx * ny);

        for (int i = 0; i < nx * ny; i++) {
            rho[i] = P_cc.rho[i];
            vx[i] = P_cc.rhovx[i] / P_cc.rho[i];
            vy[i] = P_cc.rhovy[i] / P_cc.rho[i];
            vz[i] = P_cc.rhovz[i] / P_cc.rho[i];
            Bx[i] = P_cc.Bx[i];
            By[i] = P_cc.By[i];
            Bz[i] = P_cc.Bz[i];
        }
    }


    ~PrimitiveVariablesCC() = default;

    void set(int i, int j, double rho_, double vx_, double vy_, double vz_, double Bx_, double By_, double Bz_)
    {
        rho[i * nx + j] = rho_;
        vx[i * nx + j] = vx_;
        vy[i * nx + j] = vy_;
        vz[i * nx + j] = vz_;
        Bx[i * nx + j] = Bx_;
        By[i * nx + j] = By_;
        Bz[i * nx + j] = Bz_;
    }

    void set(const ReconstructedValues& rv, int index){
        rho[index] = rv.rho;
        vx[index] = rv.vx;
        vy[index] = rv.vy;
        vz[index] = rv.vz;
        Bx[index] = rv.Bx;
        By[index] = rv.By;
        Bz[index] = rv.Bz;
    }

    void init(std::vector<double> rho_, std::vector<double> vx_, std::vector<double> vy_, std::vector<double> vz_, std::vector<double> Bx_, std::vector<double> By_, std::vector<double> Bz_){
        rho = rho_;
        vx = vx_;
        vy = vy_;
        vz = vz_;
        Bx = Bx_;
        By = By_;
        Bz = Bz_;
    }

    ReconstructedValues operator[](int index) const {
        if (index < 0 || index >= nx * ny)
            throw("Index out of range");

        return {rho[index], vx[index], vy[index], vz[index], Bx[index], By[index], Bz[index], P};
    }
};

#endif // PRIMITIVE_VARIABLES_CC_HPP_