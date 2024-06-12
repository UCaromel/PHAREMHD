#ifndef CHECK_DIV_B_HPP_
#define CHECK_DIV_B_HPP_

#include <algorithm>
#include <cmath>

#include "ConservativeVariablesCC.hpp"
//#include "WrittingUtils.hpp"

inline double CheckDivB(const ConservativeVariablesCC& Cn, double Dx, double Dy, int nghost){
    std::vector<double> DivB;

    for (int j = nghost; j < Cn.ny - nghost; ++j)
    {
        for (int i = nghost; i < Cn.nx - nghost; ++i)
        {
            DivB.push_back((Cn.Bx[j][i+1] - Cn.Bx[j][i-1])/(2*Dx) + (Cn.By[j+1][i] - Cn.By[j-1][i])/(2*Dy));
        }
    }

    //saveVectorToFile(DivB, "DivB.txt");

    auto maxAbsIter = std::max_element(DivB.begin(), DivB.end(),
                                       [](double a, double b) {
                                           return std::abs(a) < std::abs(b);
                                       });

    double maxAbsDivB = std::abs(*maxAbsIter);
    return maxAbsDivB;
}

#endif //CHECK_DIV_B_HPP_