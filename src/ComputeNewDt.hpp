#ifndef COMPUTE_NEW_DT_HPP_
#define COMPUTE_NEW_DT_HPP_

#include "ConservativeVariablesCC.hpp"
#include "EquationOfState.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>

const double sigmaCFL = 0.8;

double ComPuteNewDt(const ConservativeVariablesCC& C_cc, double Dx, double Dy, int nghost){
    std::vector<double> sum;

    for (int j = nghost; j < C_cc.ny - nghost; ++j) {
        for (int i = nghost; i < C_cc.nx - nghost; ++i) {
            double c0 = std::sqrt((gam*EosP(C_cc(i,j)))/C_cc(i,j).rho); // Sound speeds
            double cax = std::sqrt((C_cc(i,j).Bx*C_cc(i,j).Bx)/(4*M_PI*C_cc(i,j).rho)); // Alfven speeds in x
            double cay = std::sqrt((C_cc(i,j).By*C_cc(i,j).By)/(4*M_PI*C_cc(i,j).rho)); // Alfven speeds in y
            double ca = std::sqrt((C_cc(i,j).Bx*C_cc(i,j).Bx + C_cc(i,j).By*C_cc(i,j).By + C_cc(i,j).Bz*C_cc(i,j).Bz)/(4*M_PI*C_cc(i,j).rho)); // Alfven speeds
            double cfastx = std::sqrt((c0*c0 + ca*ca)*0.5 + (std::sqrt((c0*c0 + ca*ca)*(c0*c0 + ca*ca) - 4*c0*c0*cax*cax))*0.5); // Fast magnetosonic speeds in x
            double cfasty = std::sqrt((c0*c0 + ca*ca)*0.5 + (std::sqrt((c0*c0 + ca*ca)*(c0*c0 + ca*ca) - 4*c0*c0*cay*cay))*0.5); // Fast magnetosonic speeds in y
            sum.push_back(cfastx/Dx + cfasty/Dy);
        }
    }
    double max = *std::max_element(sum.begin(), sum.end());
    return sigmaCFL/max;
}


#endif //COMPUTE_NEW_DT_HPP_