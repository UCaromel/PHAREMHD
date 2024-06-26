#ifndef CONSTRAINED_TRANSPORT_HPP_
#define CONSTRAINED_TRANSPORT_HPP_

#include <vector>
#include <utility>
#include <cmath>

#include "Enums.hpp"
#include "ModularityUtils.hpp"
#include "ConservativeVariables.hpp"
#include "Interface.hpp"

std::vector<std::vector<double>> ConstrainedTransportAverage(const ConservativeVariables &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs);

std::vector<std::vector<double>> ConstrainedTransportContact(const ConservativeVariables &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs);

std::vector<std::vector<double>> UCTHLL(const ConservativeVariables &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs);

void ApplyConstrainedTransport(ConservativeVariables& Cn1, const ConservativeVariables& Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct);

#endif //CONSTRAINED_TRANSPORT_HPP_
