#include "TimeIntegrator.hpp"

#include "WrittingUtils.hpp"

ConservativeVariables EulerAdvance(const ConservativeVariables& Un, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs) {
    ConservativeVariables Un1(Un.nx, Un.ny);
    PrimitiveVariables P_cc(Un);
    ConservativeVariables FluxDifx = ComputeFluxDifferenceX(P_cc, nghost, rec, sl, rs);
    ConservativeVariables FluxDify = ComputeFluxDifferenceY(P_cc, nghost, rec, sl, rs);
    Un1 = Un - FluxDifx*(Dt/Dx) - FluxDify*(Dt/Dy);
    return Un1;
}

ConservativeVariables Euler(ConservativeVariables& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct) {
    ConservativeVariables Un1(Un.nx, Un.ny);

    UpdateGhostCells(Un, nghost, bc);

    Un1 = EulerAdvance(Un, Dx, Dy, Dt, nghost, rec, sl, rs);

    ApplyConstrainedTransport(Un1, Un, Dx, Dy, Dt, nghost, rec, sl, rs, ct); // If Un1 not needed in CT scheme (else update ghost cells before)
    
    Un1.ReconstructCenteredB(nghost);

    return Un1;
}

ConservativeVariables TVDRK2(ConservativeVariables& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct) {
    ConservativeVariables Un1(Un.nx, Un.ny);

    ConservativeVariables U1 = Euler(Un, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct);

    Un1 = Un*0.5 + Euler(U1, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct)*0.5;

    return Un1;
}

ConservativeVariables TVDRK3(ConservativeVariables& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct) {
    ConservativeVariables Un1(Un.nx, Un.ny);

    ConservativeVariables U1 = Euler(Un, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct);

    ConservativeVariables U2 = Un*0.75 + Euler(U1, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct)*0.25;

    Un1 = Un*(1.0/3.0) + Euler(U2, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct)*(2.0/3.0);

    return Un1;
}
