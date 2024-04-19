#ifndef PRIMITIVE_VARIABLES_CC_HPP_
#define PRIMITIVE_VARIABLES_CC_HPP_

#include <vector>
#include <string>
#include <cmath>

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
    const double P = (5.0/3.0)/4*M_PI;

    PrimitiveVariablesCC(int nx_, int ny_);
    PrimitiveVariablesCC(const ConservativeVariablesCC& C_cc);
    ~PrimitiveVariablesCC();

    void set(const ReconstructedValues& rv, int i, int j);
    void init(const std::vector<std::vector<double>>& rho_, const std::vector<std::vector<double>>& vx_, const std::vector<std::vector<double>>& vy_, const std::vector<std::vector<double>>& vz_, const std::vector<std::vector<double>>& Bx_, const std::vector<std::vector<double>>& By_, const std::vector<std::vector<double>>& Bz_, const double P_);
    ReconstructedValues operator()(int i, int j) const;
};

#endif // PRIMITIVE_VARIABLES_CC_HPP_
