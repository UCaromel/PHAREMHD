#include "TimeIntegrator.hpp"

ConservativeVariablesCC EulerAdvance(const ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs) {
    ConservativeVariablesCC Un1(Un.nx, Un.ny);
    PrimitiveVariablesCC P_cc(Un);
    ConservativeVariablesCC FluxDifx = ComputeFluxDifferenceX(P_cc, nghost, rec, sl, rs);
    ConservativeVariablesCC FluxDify = ComputeFluxDifferenceY(P_cc, nghost, rec, sl, rs);
    Un1 = Un - FluxDifx*(Dt/Dx) - FluxDify*(Dt/Dy);
    return Un1;
}

ConservativeVariablesCC Euler(ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct) {
    ConservativeVariablesCC Un1(Un.nx, Un.ny);

    UpdateGhostCells(Un, nghost, bc);
    Un1 = EulerAdvance(Un, Dx, Dy, Dt, nghost, rec, sl, rs);
    ApplyConstrainedTransport(Un1, Un, Dx, Dy, Dt, nghost, rec, sl, ct); // If Un1 not needed in CT scheme (else update ghost cells)

    return Un1;
}

ConservativeVariablesCC TVDRK2(ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct) {
    ConservativeVariablesCC Un1(Un.nx, Un.ny);

    ConservativeVariablesCC U1 = Euler(Un, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct);

    Un1 = Un*0.5 + Euler(U1, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct)*0.5;

    return Un1;
}

ConservativeVariablesCC TVDRK3(ConservativeVariablesCC& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct) {
    ConservativeVariablesCC Un1(Un.nx, Un.ny);

    ConservativeVariablesCC U1 = Euler(Un, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct);

    ConservativeVariablesCC U2 = Un*0.75 + Euler(U1, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct)*0.25;

    Un1 = Un*(1.0/3.0) + Euler(U2, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct)*(2.0/3.0);

    return Un1;
}
