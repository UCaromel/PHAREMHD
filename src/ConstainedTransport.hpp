#ifndef CONSTRAINED_TRANSPORT_HPP_
#define CONSTRAINED_TRANSPORT_HPP_

#include <vector>

#include <iostream>

#include "ConservativeVariablesCC.hpp"

class ConstrainedTransport
{

public:
    // Edge-centered
    std::vector<std::vector<double>> vx;
    std::vector<std::vector<double>> vy;
    std::vector<std::vector<double>> Bx;
    std::vector<std::vector<double>> By;
    std::vector<std::vector<double>> Ez;
    // Face centered
    std::vector<std::vector<double>> bx;
    std::vector<std::vector<double>> by;


    std::vector<std::vector<double>> BX;
    std::vector<std::vector<double>> BY;
    // (i,j) cell index
    ConstrainedTransport(const ConservativeVariablesCC &Cn /* Assuming ghost cells are added */, double Dx, double Dy, double Dt, int nghost)
    {

        vx.resize(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
        vy.resize(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
        Bx.resize(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
        By.resize(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
        Ez.resize(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
        bx.resize(Cn.ny - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
        by.resize(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx - 2 * nghost));

        BX.resize(Cn.ny - 2 * nghost, std::vector<double>(Cn.nx - 2 * nghost));
        BY.resize(Cn.ny - 2 * nghost, std::vector<double>(Cn.nx - 2 * nghost));

        for (int j = nghost; j <= Cn.ny - nghost; ++j)
        {
            for (int i = nghost; i <= Cn.nx - nghost; ++i)
            {
                vx[j - nghost][i - nghost] = 0.25 * ((Cn.rhovx[j][i] / Cn.rho[j][i]) + (Cn.rhovx[j][i - 1] / Cn.rho[j][i - 1]) + (Cn.rhovx[j - 1][i] / Cn.rho[j - 1][i]) + (Cn.rhovx[j - 1][i - 1] / Cn.rho[j - 1][i - 1]));
                vy[j - nghost][i - nghost] = 0.25 * ((Cn.rhovy[j][i] / Cn.rho[j][i]) + (Cn.rhovy[j][i - 1] / Cn.rho[j][i - 1]) + (Cn.rhovy[j - 1][i] / Cn.rho[j - 1][i]) + (Cn.rhovy[j - 1][i - 1] / Cn.rho[j - 1][i - 1]));
                Bx[j - nghost][i - nghost] = 0.25 * (Cn.Bx[j][i] + Cn.Bx[j][i - 1] + Cn.Bx[j - 1][i] + Cn.Bx[j - 1][i - 1]);
                By[j - nghost][i - nghost] = 0.25 * (Cn.By[j][i] + Cn.By[j][i - 1] + Cn.By[j - 1][i] + Cn.By[j - 1][i - 1]);
                Ez[j - nghost][i - nghost] = vy[j - nghost][i - nghost] * Bx[j - nghost][i - nghost] - vx[j - nghost][i - nghost] * By[j - nghost][i - nghost];
            }
        }

        for (int j = nghost; j < Cn.ny - nghost; ++j)
        {
            for (int i = nghost; i <= Cn.nx - nghost; ++i)
            {
                bx[j - nghost][i - nghost] = 0.5 * (Cn.Bx[j][i] + Cn.Bx[j][i - 1]);
                bx[j - nghost][i - nghost] -= (Dt / Dy) * (Ez[j + 1 - nghost][i - nghost] - Ez[j - nghost][i - nghost]);
            }
        }

        for (int j = nghost; j <= Cn.ny - nghost; ++j)
        {
            for (int i = nghost; i < Cn.nx - nghost; ++i)
            {
                by[j - nghost][i - nghost] = 0.5 * (Cn.By[j][i] + Cn.By[j - 1][i]);
                by[j - nghost][i - nghost] += (Dt / Dx) * (Ez[j - nghost][i + 1 - nghost] - Ez[j - nghost][i - nghost]);
            }
        }

        for (int j = 0; j < Cn.ny - 2 * nghost; ++j)
        {
            for (int i = 0; i < Cn.nx - 2 * nghost; ++i)
            {
                BX[j][i] = 0.5 * (bx[j][i + 1] + bx[j][i]);
                BY[j][i] = 0.5 * (by[j + 1][i] + by[j][i]);
            }
        }

        //std::cout << "Numerical Div : " << ((bx[10][11] - bx[10][10]) / Dx + (by[11][10] - by[10][10]) / Dy) << std::endl;
    }
    ~ConstrainedTransport() = default;
};

inline void ApplyConstrainedTransport(ConservativeVariablesCC& Cn1, const ConservativeVariablesCC& Cn, double Dx, double Dy, double Dt, int nghost){
    ConstrainedTransport CT(Cn, Dx, Dy, Dt, nghost);
    for(int j = nghost; j < Cn1.ny - nghost; ++j){
        for(int i = nghost; i < Cn1.nx - nghost; ++i){
            Cn1.Bx[j][i] = CT.BX[j - nghost][i - nghost];
            Cn1.By[j][i] = CT.BY[j - nghost][i - nghost];
        }
    }
}

#endif //CONSTRAINED_TRANSPORT_HPP_