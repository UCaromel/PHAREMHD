#ifndef CONSERVATIVE_VARIABLES_CC_HPP_
#define CONSERVATIVE_VARIABLES_CC_HPP_

#include <vector>
#include <string>

#include "PrimitiveVariables.hpp"

class PrimitiveVariables;

class ConservativeVariables {
public:
    int nx;
    int ny;
    
    // Cell-centered
    std::vector<std::vector<double>> rho;
    std::vector<std::vector<double>> rhovx;
    std::vector<std::vector<double>> rhovy;
    std::vector<std::vector<double>> rhovz;
    std::vector<std::vector<double>> Bx;
    std::vector<std::vector<double>> By;
    std::vector<std::vector<double>> Bz;
    std::vector<std::vector<double>> Etot;

    // Face-centered
    std::vector<std::vector<double>> Bxf;
    std::vector<std::vector<double>> Byf;

    ConservativeVariables(int nx, int ny);
    ConservativeVariables(const PrimitiveVariables& P_cc);
    ~ConservativeVariables();

    void ReconstructCenteredB(int nghost);
    void set(const ReconstructedValues& rv, int i, int j);
    void setflux(const ReconstructedValues& rv, int i, int j);
    ReconstructedValues operator()(int i, int j) const;
    ConservativeVariables operator*(double scalar) const;
    ConservativeVariables operator-(const ConservativeVariables& other) const;
    ConservativeVariables operator+(const ConservativeVariables& other) const;
    ConservativeVariables& operator=(const ConservativeVariables& other);
};

#endif // CONSERVATIVE_VARIABLES_CC_HPP_
