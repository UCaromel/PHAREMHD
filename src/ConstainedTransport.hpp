#ifndef CONSTRAINED_TRANSPORT_HPP_
#define CONSTRAINED_TRANSPORT_HPP_

#include <vector>

#include "ConservativeVariablesCC.hpp"

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
    ConstrainedTransport(const ConservativeVariablesCC& C_cc /* Assuming ghost cells are added */, int i /* (0 to (nx-1)) + nghost */, int j /* (0 to (ny-1)) + nghost */, double Dx, double Dy, double Dt, int nghost){

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
        by.push_back(0.5*(C_cc.By[j][i] + C_cc.By[j-1][i]));

        // Top face
        by.push_back(0.5*(C_cc.By[j][i] + C_cc.By[j+1][i]));


        double bx0 = bx[0] - (Dt/Dy)*(Ez[1]-Ez[0]);
        double bx1 = bx[1] - (Dt/Dy)*(Ez[3]-Ez[2]);
        double by0 = by[0] - (Dt/Dx)*(Ez[2]-Ez[0]);
        double by1 = by[1] - (Dt/Dx)*(Ez[3]-Ez[1]);

        BX = 0.5*(bx0 + bx1);
        BY = 0.5*(by0 + by1);
    }
    ~ConstrainedTransport() = default;
};

void ApplyConstrainedTransport(ConservativeVariablesCC& Cn1, const ConservativeVariablesCC& Cn, double Dx, double Dy, double Dt, int nghost){
    for(int j = nghost; j < Cn1.ny - nghost; j++){
        for(int i = nghost; i < Cn1.nx - nghost; i++){
            ConstrainedTransport CT(Cn, i, j, Dx, Dy, Dt, nghost);
            Cn1.Bx[j][i] = CT.BX;
            Cn1.By[j][i] = CT.BY;
        }
    }
}

#endif //CONSTRAINED_TRANSPORT_HPP_