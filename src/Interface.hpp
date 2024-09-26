#ifndef INTERFACE_HPP_
#define INTERFACE_HPP_

#include <cmath>
#include <stdexcept>
#include <vector>

#include "Enums.hpp"
#include "EquationOfState.hpp"
#include "ModularityUtils.hpp"
#include "PrimitiveVariables.hpp"
#include "ReconstructedValues.hpp"
#include "SlopeLimiter.hpp"
#include "Utils.hpp"

enum struct Dir { X, Y };

ReconstructedValues ComputeFluxVector(const ReconstructedValues &u, Dir dir);

void AddNonIdealFlux(ReconstructedValues &f, const ReconstructedValues &u,
                     double Jx, double Jy, double Jz, double LaplJx,
                     double LaplJy, double LaplJz, OptionalPhysics OptP,
                     Dir dir);
class Interface {
public:
  ReconstructedValues uL, uR, fL, fR;
  double SL, SR, Splus;
  double SLb, SRb, Splusb;
  double cfastxL, cfastxR, cfastyL, cfastyR;
  double cwxL, cwxR, cwyL, cwyR;
  OptionalPhysics OP;

  // For vector initialisation in UCT
  Interface();

  // (i,j) interface index.
  Interface(const PrimitiveVariables &P_cc, int i, int j, double Dx, double Dy,
            int nghost, Reconstruction rec, Slope sl, OptionalPhysics OptP,
            Dir dir);
  ~Interface();
};

#endif // INTERFACE_HPP_
