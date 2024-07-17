#ifndef CONSTRAINED_TRANSPORT_HPP_
#define CONSTRAINED_TRANSPORT_HPP_

#include <vector>
#include <utility>
#include <cmath>

#include "Enums.hpp"
#include "ModularityUtils.hpp"
#include "ConservativeVariables.hpp"
#include "Interface.hpp"
#include "Utils.hpp"

std::vector<std::vector<double>> ConstrainedTransportAverage(const ConservativeVariables &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs, OptionalPhysics OptP);

std::vector<std::vector<double>> ConstrainedTransportContact(const ConservativeVariables &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs, OptionalPhysics OptP);

std::vector<std::vector<double>> UCTHLL(const ConservativeVariables &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs, OptionalPhysics OptP);

void ApplyConstrainedTransport(ConservativeVariables& Cn1, const ConservativeVariables& Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct, OptionalPhysics OptP);

#endif //CONSTRAINED_TRANSPORT_HPP_
