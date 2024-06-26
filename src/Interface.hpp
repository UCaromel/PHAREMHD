#ifndef INTERFACE_HPP_
#define INTERFACE_HPP_

#include <cmath>
#include <stdexcept>

#include "ReconstructedValues.hpp"
#include "PrimitiveVariables.hpp"
#include "SlopeLimiter.hpp"
#include "EquationOfState.hpp"
#include "Enums.hpp"
#include "ModularityUtils.hpp"

enum struct Dir {
    X,
    Y
};

ReconstructedValues ComputeFluxVector(ReconstructedValues u, Dir dir);

class Interface {
public:
    ReconstructedValues uL, uR, fL, fR;
    double SL, SR, Splus;
    double cfastxL, cfastxR, cfastyL, cfastyR;

    // For vector initialisation in UCT
    Interface();

    // (i,j) interface index.
    Interface(const PrimitiveVariables& P_cc, int i, int j, Reconstruction rec, Slope sl, int nghost, Dir dir);
    ~Interface();
};

#endif // INTERFACE_HPP_