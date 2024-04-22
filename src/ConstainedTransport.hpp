#ifndef CONSTRAINED_TRANSPORT_HPP_
#define CONSTRAINED_TRANSPORT_HPP_

#include <vector>

#include "ConservativeVariablesCC.hpp"
#include "AddGhostCells.hpp"

class ConstrainedTransport
{
    // Edge-centered
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> Bx;
    std::vector<double> By;
    std::vector<double> Ez;
    // Face centered
    std::vector<double> bx;
    std::vector<double> by;

    void calculateEdge(const ConservativeVariablesCC& C_cc, int i, int j, int i2, int j2, int index) {
        vx.push_back(0.25*((C_cc.rhovx[j][i]/C_cc.rho[j][i]) + (C_cc.rhovx[j][i2]/C_cc.rho[j][i2]) 
                            + (C_cc.rhovx[j2][i]/C_cc.rho[j2][i]) + (C_cc.rhovx[j2][i2]/C_cc.rho[j2][i2])));
        vy.push_back(0.25*((C_cc.rhovy[j][i]/C_cc.rho[j][i]) + (C_cc.rhovy[j][i2]/C_cc.rho[j][i2]) 
                            + (C_cc.rhovy[j2][i]/C_cc.rho[j2][i]) + (C_cc.rhovy[j2][i2]/C_cc.rho[j2][i2])));
        Bx.push_back(0.25*(C_cc.Bx[j][i] + C_cc.Bx[j][i2] + C_cc.Bx[j2][i] + C_cc.Bx[j2][i2]));
        By.push_back(0.25*(C_cc.By[j][i] + C_cc.By[j][i2] + C_cc.By[j2][i] + C_cc.By[j2][i2]));
        Ez.push_back(vy[index]*Bx[index] - vx[index]*By[index]);
    }

public:
    double BX, BY;
    // (i,j) cell index
    ConstrainedTransport(const ConservativeVariablesCC& C_cc /* Assuming ghost cells are added */, int i /* 0 to (nx-1) */, int j /* 0 to (ny-1) */, double Dx, double Dy, double Dt, int nghost){
        i += nghost;
        j += nghost;

        // Bottom-left edge
        calculateEdge(C_cc, i, j, i-1, j-1, 0);

        // Top-left edge
        calculateEdge(C_cc, i, j, i-1, j+1, 1);

        // Bottom-right edge
        calculateEdge(C_cc, i, j, i+1, j-1, 2);

        // Top-right edge
        calculateEdge(C_cc, i, j, i+1, j+1, 3);

        // Left face
        bx.push_back(0.5*(C_cc.Bx[j][i] + C_cc.Bx[j][i-1]));

        // Right face
        bx.push_back(0.5*(C_cc.Bx[j][i] + C_cc.Bx[j][i+1]));

        // Bottom face
        by.push_back(0.5*(C_cc.Bx[j][i] + C_cc.Bx[j-1][i]));

        // Top face
        by.push_back(0.5*(C_cc.Bx[j][i] + C_cc.Bx[j+1][i]));


        bx[0] = bx[0] - (Dt/Dy)*(Ez[1]-Ez[0]);
        bx[1] = bx[1] - (Dt/Dy)*(Ez[3]-Ez[2]);
        by[0] = by[0] - (Dt/Dx)*(Ez[2]-Ez[0]);
        by[1] = by[1] - (Dt/Dx)*(Ez[3]-Ez[1]);

        BX = 0.5*(bx[0] + bx[1]);
        BY = 0.5*(by[0] + by[1]);
    }
    ~ConstrainedTransport() = default;
};

void ApplyConstrainedTransport(ConservativeVariablesCC& C_cc, double Dx, double Dy, double Dt, int nghost){
    ConservativeVariablesCC Cghost = AddGhostCells(C_cc, nghost);
    for(int j = 0; j < C_cc.ny; j++){
        for(int i = 0; i < C_cc.nx; i++){
            ConstrainedTransport CT(Cghost, i, j, Dx, Dy, Dt, nghost);
            C_cc.Bx[j][i] = CT.BX;
            C_cc.By[j][i] = CT.BY;
        }
    }
}

#endif //CONSTRAINED_TRANSPORT_HPP_