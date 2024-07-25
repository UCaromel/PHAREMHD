#include "TimeIntegrator.hpp"

#include "WrittingUtils.hpp"

ConservativeVariables EulerAdvance(const ConservativeVariables& Un, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs, OptionalPhysics OptP) {
    ConservativeVariables Un1(Un.nx, Un.ny);
    PrimitiveVariables P_cc(Un);
/*
    static int j = 0;
    std::ostringstream filenamex;
    filenamex << "whislerwaveres/" <<"Jx2_"<<j<< ".txt";
    saveVectorToFile(Un.Jx, filenamex.str());
    std::ostringstream filenamey;
    filenamey << "whislerwaveres/" <<"Jy2_"<<j<< ".txt";
    saveVectorToFile(Un.Jy, filenamey.str());
    std::ostringstream filenamez;
    filenamez << "whislerwaveres/" <<"Jz2_"<<j<< ".txt";
    saveVectorToFile(Un.Jz, filenamez.str());
    j++;
*/
    ConservativeVariables FluxDifx = ComputeFluxDifferenceX(P_cc, Dx, Dy, nghost, rec, sl, rs, OptP);
    ConservativeVariables FluxDify = ComputeFluxDifferenceY(P_cc, Dx, Dy, nghost, rec, sl, rs, OptP);
    Un1 = Un - FluxDifx*(Dt/Dx) - FluxDify*(Dt/Dy);
    return Un1;
}

ConservativeVariables Euler(ConservativeVariables& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct, OptionalPhysics OptP) {
    ConservativeVariables Un1(Un.nx, Un.ny);

    UpdateGhostCells(Un, nghost, bc);

    if (OptP == OptionalPhysics::HallResHyper) {
        ComputeJ(Un, Dx, Dy, nghost);
        UpdateGhostJ(Un, nghost, bc);
        static int j = 0;
        std::ostringstream filenamex;
        filenamex << "whislerwaveres/" <<"Jx_"<<j<< ".txt";
        saveVectorToFile(Un.Jx, filenamex.str());
        std::ostringstream filenamey;
        filenamey << "whislerwaveres/" <<"Jy_"<<j<< ".txt";
        saveVectorToFile(Un.Jy, filenamey.str());
        std::ostringstream filenamez;
        filenamez << "whislerwaveres/" <<"Jz_"<<j<< ".txt";
        saveVectorToFile(Un.Jz, filenamez.str());
        j++;
    }

    Un1 = EulerAdvance(Un, Dx, Dy, Dt, nghost, rec, sl, rs, OptP);

    //ApplyConstrainedTransport(Un1, Un, Dx, Dy, Dt, nghost, rec, sl, rs, ct, OptP); // If Un1 not needed in CT scheme (else update ghost cells before)
    
    //Un1.ReconstructCenteredB(nghost);

    return Un1;
}

ConservativeVariables TVDRK2(ConservativeVariables& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct, OptionalPhysics OptP) {
    ConservativeVariables Un1(Un.nx, Un.ny);

    ConservativeVariables U1 = Euler(Un, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct, OptP);

    Un1 = Un*0.5 + Euler(U1, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct, OptP)*0.5;

    return Un1;
}

ConservativeVariables TVDRK3(ConservativeVariables& Un, double Dx, double Dy, double Dt, int nghost, BoundaryConditions bc, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct, OptionalPhysics OptP) {
    ConservativeVariables Un1(Un.nx, Un.ny);

    ConservativeVariables U1 = Euler(Un, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct, OptP);

    ConservativeVariables U2 = Un*0.75 + Euler(U1, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct, OptP)*0.25;

    Un1 = Un*(1.0/3.0) + Euler(U2, Dx, Dy, Dt, nghost, bc, rec, sl, rs, ct, OptP)*(2.0/3.0);

    return Un1;
}
