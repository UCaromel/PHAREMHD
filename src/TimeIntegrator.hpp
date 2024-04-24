#ifndef TIME_INTEGRATOR_HPP_
#define TIME_INTEGRATOR_HPP_

#include <vector>

#include "ConservativeVariablesCC.hpp"
#include "GodunovFlux.hpp"
#include "ConstainedTransport.hpp"
#include "AddGhostCells.hpp"

ConservativeVariablesCC EulerAdvance(const ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int order, int nghost){
    ConservativeVariablesCC Un1(Un.nx, Un.ny);
    PrimitiveVariablesCC P_cc(Un);
    ConservativeVariablesCC FluxDifx = ComputeFluxDifferenceX(P_cc, order, nghost);
    ConservativeVariablesCC FluxDify = ComputeFluxDifferenceY(P_cc, order, nghost);
    Un1 = Un - FluxDifx*(Dt/Dx) - FluxDify*(Dt/Dy);
    return Un1;
}

ConservativeVariablesCC Euler(ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int order, int nghost){
    ConservativeVariablesCC Un1(Un.nx, Un.ny);

    UpdateGhostCells(Un, nghost);
    Un1 = EulerAdvance(Un, Dx, Dy, Dt, order, nghost);
    ApplyConstrainedTransport(Un1, Un, Dx, Dy, Dt, nghost); // If Un1 not needed in CT scheme (else update ghost cells)

    return Un1;
}

ConservativeVariablesCC TVDRK2(ConservativeVariablesCC& Un,  double Dx, double Dy, double Dt, int order, int nghost){
    ConservativeVariablesCC Un1(Un.nx, Un.ny);

    UpdateGhostCells(Un, nghost);
    ConservativeVariablesCC U1 = EulerAdvance(Un, Dx, Dy, Dt, order, nghost);
    ApplyConstrainedTransport(U1, Un, Dx, Dy, Dt, nghost);

    UpdateGhostCells(U1, nghost);
    Un1 = Un*0.5 + EulerAdvance(U1, Dx, Dy, Dt, order, nghost)*0.5;
    ApplyConstrainedTransport(Un1, U1, Dx, Dy, Dt, nghost);

    return Un1;
}

ConservativeVariablesCC TVDRK3(ConservativeVariablesCC& Un,  double Dx, double Dy, double Dt, int order, int nghost){
    ConservativeVariablesCC Un1(Un.nx, Un.ny);

    UpdateGhostCells(Un, nghost);
    ConservativeVariablesCC U1 = EulerAdvance(Un, Dx, Dy, Dt, order, nghost);
    ApplyConstrainedTransport(U1, Un, Dx, Dy, Dt, nghost);

    UpdateGhostCells(U1, nghost);
    ConservativeVariablesCC U2 = Un*0.75 + EulerAdvance(U1, Dx, Dy, Dt, order, nghost)*0.25;
    ApplyConstrainedTransport(U2, U1, Dx, Dy, Dt, nghost);

    UpdateGhostCells(U2, nghost);
    Un1 = Un*(1.0/3.0) + EulerAdvance(U2, Dx, Dy, Dt, order, nghost)*(2.0/3.0);
    ApplyConstrainedTransport(Un1, U2, Dx, Dy, Dt, nghost);

    return Un1;
}


#endif //TIME_INTEGRATOR_HPP_