#ifndef CHECK_DIV_B_HPP_
#define CHECK_DIV_B_HPP_

#include <algorithm>
#include <cmath>

#include "ConservativeVariables.hpp"

inline double CheckDivB(const ConservativeVariables& Cn, double Dx, double Dy, int nghost){
    std::vector<double> DivB;

    for (int j = nghost; j < Cn.ny - nghost; ++j)
    {
        for (int i = nghost; i < Cn.nx - nghost; ++i)
        {
            DivB.push_back((Cn.Bxf[j][i+1] - Cn.Bxf[j][i])/(Dx) + (Cn.Byf[j+1][i] - Cn.Byf[j][i])/(Dy));
        }
    }

    auto maxAbsIter = std::max_element(DivB.begin(), DivB.end(),
                                       [](double a, double b) {
                                           return std::abs(a) < std::abs(b);
                                       });

    double maxAbsDivB = std::abs(*maxAbsIter);
    return maxAbsDivB;
}

#endif //CHECK_DIV_B_HPP_