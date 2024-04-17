#ifndef TIME_INTEGRATOR_HPP_
#define TIME_INTEGRATOR_HPP_

#include <vector>

#include "ConservativeVariablesCC.hpp"
#include "ReconstructedValues.hpp"
#include "GodunovFlux.hpp"

ConservativeVariablesCC EulerAdvance(const ConservativeVariablesCC& Un, ConservativeVariablesCC& FluxDif, double Dx, double Dt){
    ConservativeVariablesCC Un1(Un.nx, Un.ny);
    Un1 = Un - FluxDif*(Dt/Dx);
}

ConservativeVariablesCC TVDRK2(const ConservativeVariablesCC& Un, ConservativeVariablesCC& FluxDif, double Dx, double Dt, int order){
    ConservativeVariablesCC Un1(Un.nx, Un.ny);
    ConservativeVariablesCC U1 = EulerAdvance(Un, FluxDif, Dx, Dt);
    PrimitiveVariablesCC P_cc(U1);
    ConservativeVariablesCC FluxDif1 = ComputeFluxDifference(P_cc, order);
    //TODO
}


#endif //TIME_INTEGRATOR_HPP_