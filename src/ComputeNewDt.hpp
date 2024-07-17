#ifndef COMPUTE_NEW_DT_HPP_
#define COMPUTE_NEW_DT_HPP_

#include "ConservativeVariables.hpp"
#include "EquationOfState.hpp"
#include "Enums.hpp"
#include "PhysicalConstants.hpp"

#include <cmath>
#include <algorithm>

inline double ComPuteNewDt(const ConservativeVariables& C_cc, double Dx, double Dy, int nghost, OptionalPhysics OptP){
    const PhysicalConstants& pc = PhysicalConstants::getInstance();
    std::vector<double> sum;

    for (int j = nghost; j < C_cc.ny - nghost; ++j) {
        for (int i = nghost; i < C_cc.nx - nghost; ++i) {

            double c0 = std::sqrt((pc.gam*EosP(C_cc(i,j)))/C_cc(i,j).rho); // Sound speeds
            double cax = std::sqrt((C_cc(i,j).Bx*C_cc(i,j).Bx)/(C_cc(i,j).rho)); // Alfven speeds in x
            double cay = std::sqrt((C_cc(i,j).By*C_cc(i,j).By)/(C_cc(i,j).rho)); // Alfven speeds in y
            double ca = std::sqrt((C_cc(i,j).Bx*C_cc(i,j).Bx + C_cc(i,j).By*C_cc(i,j).By + C_cc(i,j).Bz*C_cc(i,j).Bz)/(C_cc(i,j).rho)); // Alfven speeds
            double cfastx = std::sqrt((c0*c0 + ca*ca)*0.5 + (std::sqrt((c0*c0 + ca*ca)*(c0*c0 + ca*ca) - 4*c0*c0*cax*cax))*0.5); // Fast magnetosonic speeds in x
            double cfasty = std::sqrt((c0*c0 + ca*ca)*0.5 + (std::sqrt((c0*c0 + ca*ca)*(c0*c0 + ca*ca) - 4*c0*c0*cay*cay))*0.5); // Fast magnetosonic speeds in y

            if (OptP == OptionalPhysics::HallResHyper) {
                double cwx = std::abs(C_cc(i,j).Bx) * M_PI / (C_cc(i,j).rho * Dx);
                double cwy = std::abs(C_cc(i,j).By) * M_PI / (C_cc(i,j).rho * Dy);

                cfastx += cwx;
                cfasty += cwy;
            }

            double dtij = (C_cc(i,j).vx/C_cc(i,j).rho + cfastx)/Dx + (C_cc(i,j).vy/C_cc(i,j).rho + cfasty)/Dy;

            if (OptP == OptionalPhysics::HallResHyper) {
                dtij += 2.0*pc.eta * (1.0/(Dx*Dx) + 1.0/(Dy*Dy));
            }
            
            sum.push_back(dtij);
        }
    }
    double max = *std::max_element(sum.begin(), sum.end());
    return pc.sigmaCFL/max;
}


#endif //COMPUTE_NEW_DT_HPP_