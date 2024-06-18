#ifndef CONSTRAINED_TRANSPORT_HPP_
#define CONSTRAINED_TRANSPORT_HPP_

#include <vector>
#include <utility>
#include <cmath>

#include "Enums.hpp"
#include "ModularityUtils.hpp"
#include "ConservativeVariablesCC.hpp"
#include "Interface.hpp"

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> ConstrainedTransportAverage(const ConservativeVariablesCC &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs);

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> ConstrainedTransportContact(const ConservativeVariablesCC &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs);

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> UCTHLL(const ConservativeVariablesCC &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs);

void ApplyConstrainedTransport(ConservativeVariablesCC& Cn1, const ConservativeVariablesCC& Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct);

#endif //CONSTRAINED_TRANSPORT_HPP_
