#ifndef TIME_INTEGRATOR_HPP_
#define TIME_INTEGRATOR_HPP_

#include <vector>

#include "ConservativeVariablesCC.hpp"
#include "GodunovFlux.hpp"
#include "ConstainedTransport.hpp"
#include "AddGhostCells.hpp"

inline ConservativeVariablesCC EulerAdvance(const ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int order, int nghost){
    ConservativeVariablesCC Un1(Un.nx, Un.ny);
    PrimitiveVariablesCC P_cc(Un);
    ConservativeVariablesCC FluxDifx = ComputeFluxDifferenceX(P_cc, order, nghost);
    ConservativeVariablesCC FluxDify = ComputeFluxDifferenceY(P_cc, order, nghost);
    Un1 = Un - FluxDifx*(Dt/Dx) - FluxDify*(Dt/Dy);
    return Un1;
}

inline ConservativeVariablesCC Euler(ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int order, int nghost){
    ConservativeVariablesCC Un1(Un.nx, Un.ny);

    UpdateGhostCells(Un, nghost);
    Un1 = EulerAdvance(Un, Dx, Dy, Dt, order, nghost);
    ApplyConstrainedTransport(Un1, Un, Dx, Dy, Dt, nghost); // If Un1 not needed in CT scheme (else update ghost cells)

    return Un1;
}

inline ConservativeVariablesCC TVDRK2(ConservativeVariablesCC& Un,  double Dx, double Dy, double Dt, int order, int nghost){
    ConservativeVariablesCC Un1(Un.nx, Un.ny);

    ConservativeVariablesCC U1 = Euler(Un, Dx, Dy, Dt, order, nghost);

    Un1 = Un*0.5 + Euler(U1, Dx, Dy, Dt, order, nghost)*0.5;

    return Un1;
}

inline ConservativeVariablesCC TVDRK3(ConservativeVariablesCC& Un,  double Dx, double Dy, double Dt, int order, int nghost){
    ConservativeVariablesCC Un1(Un.nx, Un.ny);

    ConservativeVariablesCC U1 = Euler(Un, Dx, Dy, Dt, order, nghost);

    ConservativeVariablesCC U2 = Un*0.75 + Euler(U1, Dx, Dy, Dt, order, nghost)*0.25;

    Un1 = Un*(1.0/3.0) + Euler(U2, Dx, Dy, Dt, order, nghost)*(2.0/3.0);

    return Un1;
}


#endif //TIME_INTEGRATOR_HPP_