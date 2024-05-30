#include "TimeIntegrator.hpp"

ConservativeVariablesCC EulerAdvance(const ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Riemann rs) {
    ConservativeVariablesCC Un1(Un.nx, Un.ny);
    PrimitiveVariablesCC P_cc(Un);
    ConservativeVariablesCC FluxDifx = ComputeFluxDifferenceX(P_cc, nghost, rec, rs);
    ConservativeVariablesCC FluxDify = ComputeFluxDifferenceY(P_cc, nghost, rec, rs);
    Un1 = Un - FluxDifx*(Dt/Dx) - FluxDify*(Dt/Dy);
    return Un1;
}

ConservativeVariablesCC Euler(ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Riemann rs, CTMethod ct) {
    ConservativeVariablesCC Un1(Un.nx, Un.ny);

    UpdateGhostCells(Un, nghost);
    Un1 = EulerAdvance(Un, Dx, Dy, Dt, nghost, rec, rs);
    ApplyConstrainedTransport(Un1, Un, Dx, Dy, Dt, nghost, rec, ct); // If Un1 not needed in CT scheme (else update ghost cells)

    return Un1;
}

ConservativeVariablesCC TVDRK2(ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Riemann rs, CTMethod ct) {
    ConservativeVariablesCC Un1(Un.nx, Un.ny);

    ConservativeVariablesCC U1 = Euler(Un, Dx, Dy, Dt, nghost, rec, rs, ct);

    Un1 = Un*0.5 + Euler(U1, Dx, Dy, Dt, nghost, rec, rs, ct)*0.5;

    return Un1;
}

ConservativeVariablesCC TVDRK3(ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Riemann rs, CTMethod ct) {
    ConservativeVariablesCC Un1(Un.nx, Un.ny);

    ConservativeVariablesCC U1 = Euler(Un, Dx, Dy, Dt, nghost, rec, rs, ct);

    ConservativeVariablesCC U2 = Un*0.75 + Euler(U1, Dx, Dy, Dt, nghost, rec, rs, ct)*0.25;

    Un1 = Un*(1.0/3.0) + Euler(U2, Dx, Dy, Dt, nghost, rec, rs, ct)*(2.0/3.0);

    return Un1;
}
