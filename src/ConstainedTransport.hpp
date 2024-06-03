#ifndef CONSTRAINED_TRANSPORT_HPP_
#define CONSTRAINED_TRANSPORT_HPP_

#include <vector>
#include <utility>
#include <cmath>

#include "ModularityUtils.hpp"
#include "ConservativeVariablesCC.hpp"

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> ConstrainedTransportAverage(const ConservativeVariablesCC &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl);

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> UCTHLL(const ConservativeVariablesCC &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl);

void ApplyConstrainedTransport(ConservativeVariablesCC& Cn1, const ConservativeVariablesCC& Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, CTMethod ct);

#endif //CONSTRAINED_TRANSPORT_HPP_
