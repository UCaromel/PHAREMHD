#ifndef INTERFACE_HPP_
#define INTERFACE_HPP_

#include <cmath>
#include <stdexcept>

#include "ReconstructedValues.hpp"
#include "PrimitiveVariablesCC.hpp"
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
    Interface(const PrimitiveVariablesCC& P_cc /* Assuming ghost cells are added */, int i /* (0 to nx) + nghost */, int j /* (0 to ny) + nghost */, Reconstruction rec, Slope sl, int nghost, Dir dir);
    ~Interface();
};

#endif // INTERFACE_HPP_