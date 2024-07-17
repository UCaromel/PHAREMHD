#ifndef INTERFACE_HPP_
#define INTERFACE_HPP_

#include <cmath>
#include <stdexcept>
#include <vector>

#include "ReconstructedValues.hpp"
#include "PrimitiveVariables.hpp"
#include "SlopeLimiter.hpp"
#include "EquationOfState.hpp"
#include "Enums.hpp"
#include "ModularityUtils.hpp"
#include "Utils.hpp"

enum struct Dir {
    X,
    Y
};

ReconstructedValues ComputeFluxVector(const ReconstructedValues& u, Dir dir);
void AddNonIdealFlux(ReconstructedValues& f, const ReconstructedValues& u, double Jz, double LaplJz, OptionalPhysics OptP, Dir dir);
std::pair<std::pair<double, double>, std::pair<double, double>> ComputeRiemannJ(const std::vector<std::vector<double>>& J, int i, int j, double Dx, double Dy, Reconstruction rec, Slope sl, Dir Dir);

class Interface {
public:
    ReconstructedValues uL, uR, fL, fR;
    double SL, SR, Splus;
    double cfastxL, cfastxR, cfastyL, cfastyR;
    double cwxL, cwxR, cwyL, cwyR;

    // For vector initialisation in UCT
    Interface();

    // (i,j) interface index.
    Interface(const PrimitiveVariables& P_cc, int i, int j, double Dx, double Dy, int nghost, Reconstruction rec, Slope sl, OptionalPhysics OptP, Dir dir);
    ~Interface();
};

#endif // INTERFACE_HPP_