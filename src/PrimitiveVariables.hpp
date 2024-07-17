#ifndef PRIMITIVE_VARIABLES_HPP_
#define PRIMITIVE_VARIABLES_HPP_

#include <vector>
#include <string>
#include <cmath>

#include "ReconstructedValues.hpp"
#include "ConservativeVariables.hpp"

class ConservativeVariables;

class PrimitiveVariables {
public:
    int nx;
    int ny;

    // Cell-centered
    std::vector<std::vector<double>> rho;
    std::vector<std::vector<double>> vx;
    std::vector<std::vector<double>> vy;
    std::vector<std::vector<double>> vz;
    std::vector<std::vector<double>> Bx;
    std::vector<std::vector<double>> By;
    std::vector<std::vector<double>> Bz;
    std::vector<std::vector<double>> P;

    // Face-centered
    std::vector<std::vector<double>> Bxf;
    std::vector<std::vector<double>> Byf;

    std::vector<std::vector<double>> Jx;
    std::vector<std::vector<double>> Jy;

    // Edge-centered
    std::vector<std::vector<double>> Jz;

    PrimitiveVariables(int nx_, int ny_);
    PrimitiveVariables(const ConservativeVariables& C_cc);
    ~PrimitiveVariables();

    void ReconstructCenteredB(int nghost);
    void set(const ReconstructedValues& rv, int i, int j);
    void init(const std::vector<std::vector<double>>& rho_, const std::vector<std::vector<double>>& vx_, const std::vector<std::vector<double>>& vy_, const std::vector<std::vector<double>>& vz_, const std::vector<std::vector<double>>& Bx_, const std::vector<std::vector<double>>& By_, const std::vector<std::vector<double>>& Bz_, const std::vector<std::vector<double>>& P_);
    ReconstructedValues operator()(int i, int j) const;
    PrimitiveVariables& operator=(const PrimitiveVariables& other);
};

#endif // PRIMITIVE_VARIABLES_HPP_
