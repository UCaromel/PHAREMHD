#ifndef EQUATION_OF_STATE_HPP_
#define EQUATION_OF_STATE_HPP_

#include "ReconstructedValues.hpp"
#include "PhysicalConstants.hpp"

double EosEtot(const ReconstructedValues& rv);
double EosP(const ReconstructedValues& rv);

#endif // EQUATION_OF_STATE_HPP_
