#ifndef CONSERVATIVE_VARIABLES_CC_HPP_
#define CONSERVATIVE_VARIABLES_CC_HPP_

#include <vector>
#include <string>

#include "PrimitiveVariablesCC.hpp"

class PrimitiveVariablesCC;

class ConservativeVariablesCC {
public:
    int nx;
    int ny;
    std::vector<std::vector<double>> rho;
    std::vector<std::vector<double>> rhovx;
    std::vector<std::vector<double>> rhovy;
    std::vector<std::vector<double>> rhovz;
    std::vector<std::vector<double>> Bx;
    std::vector<std::vector<double>> By;
    std::vector<std::vector<double>> Bz;
    const double P = (5.0/3.0)/4*M_PI;

    ConservativeVariablesCC(int nx, int ny);
    ConservativeVariablesCC(const PrimitiveVariablesCC& P_cc);
    ~ConservativeVariablesCC();
    
    void set(const ReconstructedValues& rv, int i, int j);
    ReconstructedValues operator()(int i, int j) const;
    ConservativeVariablesCC operator*(double scalar) const;
    ConservativeVariablesCC operator-(const ConservativeVariablesCC& other) const;
    ConservativeVariablesCC operator+(const ConservativeVariablesCC& other) const;
    ConservativeVariablesCC& operator=(const ConservativeVariablesCC& other);
};

#endif // CONSERVATIVE_VARIABLES_CC_HPP_
