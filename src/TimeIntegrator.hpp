#ifndef TIME_INTEGRATOR_HPP_
#define TIME_INTEGRATOR_HPP_

#include <vector>

#include "ConservativeVariablesCC.hpp"
#include "GodunovFlux.hpp"

ConservativeVariablesCC EulerAdvance(const ConservativeVariablesCC& Un, double Dx, double Dt, int order){
    ConservativeVariablesCC Un1(Un.nx, Un.ny);
    PrimitiveVariablesCC P_cc(Un);
    ConservativeVariablesCC FluxDif = ComputeFluxDifference(P_cc, order);
    Un1 = Un - FluxDif*(Dt/Dx);
    // Apply CT
}

ConservativeVariablesCC TVDRK2(const ConservativeVariablesCC& Un,  double Dx, double Dt, int order){
    ConservativeVariablesCC Un1(Un.nx, Un.ny);
    ConservativeVariablesCC U1 = EulerAdvance(Un, Dx, Dt, order);
    PrimitiveVariablesCC P_cc(U1);
    ConservativeVariablesCC FluxDif1 = ComputeFluxDifference(P_cc, order);
    Un1 = Un*0.5 + EulerAdvance(U1, Dx, Dt, order)*0.5;
}


#endif //TIME_INTEGRATOR_HPP_