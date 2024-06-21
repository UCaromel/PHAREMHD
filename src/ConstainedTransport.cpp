#include "ConstainedTransport.hpp"

std::vector<std::vector<double>> ConstrainedTransportAverage(const ConservativeVariables &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs)
{
    // Edge-centered
    std::vector<std::vector<double>> vx(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
    std::vector<std::vector<double>> vy(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
    std::vector<std::vector<double>> Bx(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
    std::vector<std::vector<double>> By(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
    std::vector<std::vector<double>> Ez(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));

    for (int j = nghost; j <= Cn.ny - nghost; ++j)
    {
        for (int i = nghost; i <= Cn.nx - nghost; ++i)
        {
            vx[j - nghost][i - nghost] = 0.25 * ((Cn.rhovx[j][i] / Cn.rho[j][i]) + (Cn.rhovx[j][i - 1] / Cn.rho[j][i - 1]) + (Cn.rhovx[j - 1][i] / Cn.rho[j - 1][i]) + (Cn.rhovx[j - 1][i - 1] / Cn.rho[j - 1][i - 1]));
            vy[j - nghost][i - nghost] = 0.25 * ((Cn.rhovy[j][i] / Cn.rho[j][i]) + (Cn.rhovy[j][i - 1] / Cn.rho[j][i - 1]) + (Cn.rhovy[j - 1][i] / Cn.rho[j - 1][i]) + (Cn.rhovy[j - 1][i - 1] / Cn.rho[j - 1][i - 1]));
            Bx[j - nghost][i - nghost] = 0.5 * (Cn.Bxf[j - 1][i] + Cn.Bxf[j][i]);
            By[j - nghost][i - nghost] = 0.5 * (Cn.Byf[j][i - 1] + Cn.Byf[j][i]);
            Ez[j - nghost][i - nghost] = vy[j - nghost][i - nghost] * Bx[j - nghost][i - nghost] - vx[j - nghost][i - nghost] * By[j - nghost][i - nghost];
        }
    }

    return Ez;
}

std::vector<std::vector<double>> ConstrainedTransportContact(const ConservativeVariables &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs){
/*    PrimitiveVariables Pn(Cn);
    std::vector<std::vector<Interface>> InterfacesX(Pn.ny + 2 - 2 * nghost, std::vector<Interface>(Pn.nx + 1 - 2 * nghost));
    std::vector<std::vector<Interface>> InterfacesY(Pn.ny + 1 - 2 * nghost, std::vector<Interface>(Pn.nx + 2 - 2 * nghost));

    std::vector<std::vector<double>> bx(Cn.ny - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
    std::vector<std::vector<double>> by(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx - 2 * nghost));

    std::vector<std::vector<double>> BX(Cn.ny - 2 * nghost, std::vector<double>(Cn.nx - 2 * nghost));
    std::vector<std::vector<double>> BY(Cn.ny - 2 * nghost, std::vector<double>(Cn.nx - 2 * nghost));

    std::vector<std::vector<double>> Ez(Pn.ny + 1 - 2 * nghost, std::vector<double>(Pn.nx + 1 - 2 * nghost));

    RiemannSolverFunction ChosenRiemannSolver = getRiemannSolver(rs);

    for (int j = nghost - 1; j < Cn.ny - nghost + 1; ++j)
    {
        for (int i = nghost; i <= Cn.nx - nghost; ++i)
        {
            InterfacesX[j - nghost + 1][i - nghost] = Interface(Pn, i, j, rec, sl, nghost, Dir::X);
        }
    }

    for (int j = nghost; j <= Cn.ny - nghost; ++j)
    {
        for (int i = nghost - 1; i < Cn.nx - nghost + 1; ++i)
        {
            InterfacesY[j - nghost][i - nghost + 1] = Interface(Pn, i, j, rec, sl, nghost, Dir::Y);
        }
    }

    for (int j = 0; j <= Cn.ny - 2*nghost; ++j)
    {
        for (int i = 0; i <= Cn.nx - 2*nghost; ++i)
        {
            Interface x = InterfacesX[j][i];
            Interface x1 = InterfacesX[j + 1][i];
            Interface y = InterfacesY[j][i];
            Interface y1 = InterfacesY[j][i + 1];

            // Used to compute sign of contact discontinuity
            double Fvx = ChosenRiemannSolver(x).rho;
            double Fvx1 = ChosenRiemannSolver(x1).rho;
            double Fvy = ChosenRiemannSolver(y).rho;
            double Fvy1 = ChosenRiemannSolver(y1).rho;

            // Rusanov fluxes for Bx and By
            double FBy_x = ChosenRiemannSolver(x).By;
            double FBy_x1 = ChosenRiemannSolver(x1).By;
            double FBx_y = ChosenRiemannSolver(y).Bx;
            double FBx_y1 = ChosenRiemannSolver(y1).Bx;
            
            // Conditionnals required by CT contact
            double sx = Fvx > 0 ? 1 : Fvx < 0 ? -1 : 0;
            double sx1 = Fvx1 > 0 ? 1 : Fvx1 < 0 ? -1 : 0;
            double sy = Fvy > 0 ? 1 : Fvy < 0 ? -1 : 0;
            double sy1 = Fvy1 > 0 ? 1 : Fvy1 < 0 ? -1 : 0;

            // Cell-centered F vector = cell-centered EMF
            int offset = nghost - 1;
            double FcBy_x = ComputeFluxVector(Pn(i + offset, j + offset), Dir::X).By;
            double FcBy_x1 = ComputeFluxVector(Pn(i + offset, j + 1 + offset), Dir::X).By;
            double FcBx_y = ComputeFluxVector(Pn(i + offset, j + offset), Dir::Y).Bx;
            double FcBx_y1 = ComputeFluxVector(Pn(i + 1 + offset, j + offset), Dir::Y).Bx;

            double FcBy_x1y1 = ComputeFluxVector(Pn(i + 1 + offset, j + 1 + offset), Dir::X).By;
            double FcBx_y1x1 = ComputeFluxVector(Pn(i + 1 + offset, j + 1 + offset), Dir::Y).Bx;

            // Compute Ez derivatives in each directions
            double dyEzS = (1.0 + sx) * (FBx_y - FcBx_y) + (1.0 - sx) * (FBx_y1 - FcBx_y1);
            double dyEzN = (1.0 + sx1) * (FBx_y - FcBx_y1) + (1.0 - sx1) * (FBx_y1 - FcBx_y1x1);
            double dxEzW = (1.0 + sy) * (FBy_x - FcBy_x) + (1.0 - sy) * (FBy_x1 - FcBy_x1);
            double dxEzE = (1.0 + sy1) * (FBy_x - FcBy_x1) + (1.0 - sy1) * (FBy_x1 - FcBy_x1y1);
           
            Ez[j][i] = 0.25*(- FBy_x - FBy_x1 + FBx_y + FBx_y1) ;//+ 0.125 * (dyEzS - dyEzN - dxEzW + dxEzE);
            if ((i == Cn.nx - 2*nghost && j == Cn.ny - 2*nghost) && !std::isnan(FBy_x)){
                std::cout<<FBy_x<<" "<<FBy_x1<<" "<<FBx_y<<" "<<FBx_y1<<" "<<Ez[j][i]<<std::endl;
            }
        }
    }

    for (int j = nghost; j < Cn.ny - nghost; ++j)
    {
        for (int i = nghost; i <= Cn.nx - nghost; ++i)
        {
            bx[j - nghost][i - nghost] = 0.5 * (Cn.Bx[j][i] + Cn.Bx[j][i - 1]);
            bx[j - nghost][i - nghost] -= (Dt / Dy) * (Ez[j + 1 - nghost][i - nghost] - Ez[j - nghost][i - nghost]);
        }
    }

    for (int j = nghost; j <= Cn.ny - nghost; ++j)
    {
        for (int i = nghost; i < Cn.nx - nghost; ++i)
        {
            by[j - nghost][i - nghost] = 0.5 * (Cn.By[j][i] + Cn.By[j - 1][i]);
            by[j - nghost][i - nghost] += (Dt / Dx) * (Ez[j - nghost][i + 1 - nghost] - Ez[j - nghost][i - nghost]);
        }
    }

    for (int j = 0; j < Cn.ny - 2 * nghost; ++j)
    {
        for (int i = 0; i < Cn.nx - 2 * nghost; ++i)
        {
            BX[j][i] = 0.5 * (bx[j][i + 1] + bx[j][i]);
            BY[j][i] = 0.5 * (by[j + 1][i] + by[j][i]);
        }
    }

    return {BX, BY};*/
}


std::vector<std::vector<double>> UCTHLL(const ConservativeVariables &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs){
    PrimitiveVariables Pn(Cn);
    std::vector<std::vector<Interface>> InterfacesX(Pn.ny + 2 - 2 * nghost, std::vector<Interface>(Pn.nx + 1 - 2 * nghost));
    std::vector<std::vector<Interface>> InterfacesY(Pn.ny + 1 - 2 * nghost, std::vector<Interface>(Pn.nx + 2 - 2 * nghost));

    std::vector<std::vector<double>> Ez(Pn.ny + 1 - 2 * nghost, std::vector<double>(Pn.nx + 1 - 2 * nghost));

    for (int j = nghost - 1; j < Cn.ny - nghost + 1; ++j)
    {
        for (int i = nghost; i <= Cn.nx - nghost; ++i)
        {
            InterfacesX[j - nghost + 1][i - nghost] = Interface(Pn, i, j, rec, sl, nghost, Dir::X);
        }
    }

    for (int j = nghost; j <= Cn.ny - nghost; ++j)
    {
        for (int i = nghost - 1; i < Cn.nx - nghost + 1; ++i)
        {
            InterfacesY[j - nghost][i - nghost + 1] = Interface(Pn, i, j, rec, sl, nghost, Dir::Y);
        }
    }

    for (int j = 0; j <= Cn.ny - 2*nghost; ++j){
        for (int i = 0; i <= Cn.nx - 2*nghost; ++i){
            Interface x = InterfacesX[j][i];
            Interface x1 = InterfacesX[j + 1][i];
            Interface y = InterfacesY[j][i];
            Interface y1 = InterfacesY[j][i + 1];

            // Compute solver's coeficients
            double alphaxR = std::max(0.0, x.uR.vx / x.uR.rho + x.cfastxR);
            double alphaxL = - std::min(0.0, x.uL.vx / x.uL.rho - x.cfastxL);

            double alphax1R = std::max(0.0, x1.uR.vx / x1.uR.rho + x1.cfastxR);
            double alphax1L = - std::min(0.0, x1.uL.vx / x1.uL.rho - x1.cfastxL);

            double alphayR = std::max(0.0, y.uR.vy / y.uR.rho + y.cfastyR);
            double alphayL = - std::min(0.0, y.uL.vy / y.uL.rho - y.cfastyL);

            double alphay1R = std::max(0.0, y1.uR.vy / y1.uR.rho + y1.cfastyR);
            double alphay1L = - std::min(0.0, y1.uL.vy / y1.uL.rho - y1.cfastyL);


            double axR = alphaxL / (alphaxR + alphaxL);
            double axL = alphaxR / (alphaxR + alphaxL);

            double ax1R = alphax1L / (alphax1R + alphax1L);
            double ax1L = alphax1R / (alphax1R + alphax1L);

            double ayR = alphayL / (alphayR + alphayL);
            double ayL = alphayR / (alphayR + alphayL);

            double ay1R = alphay1L / (alphay1R + alphay1L);
            double ay1L = alphay1R / (alphay1R + alphay1L);


            double dxR = (alphaxR*alphaxL) / (alphaxR + alphaxL);
            double dxL = (alphaxR*alphaxL) / (alphaxR + alphaxL);

            double dx1R = (alphax1R*alphax1L) / (alphax1R + alphax1L);
            double dx1L = (alphax1R*alphax1L) / (alphax1R + alphax1L);

            double dyR = (alphayR*alphayL) / (alphayR + alphayL);
            double dyL = (alphayR*alphayL) / (alphayR + alphayL);

            double dy1R = (alphay1R*alphay1L) / (alphay1R + alphay1L);
            double dy1L = (alphay1R*alphay1L) / (alphay1R + alphay1L);

            // Averaging
            double aW = 0.5 * (axL + ax1L);
            double aE = 0.5 * (axR + ax1R);
            double aS = 0.5 * (ayL + ay1L);
            double aN = 0.5 * (ayR + ay1R);

            double dW = 0.5 * (dxL + dx1L);
            double dE = 0.5 * (dxR + dx1R);
            double dS = 0.5 * (dyL + dy1L);
            double dN = 0.5 * (dyR + dy1R);

            double vxW = 0.5 * (y.uL.vx / y.uL.rho + y.uR.vx / y.uR.rho);
            double vxE = 0.5 * (y1.uL.vx / y1.uL.rho + y1.uR.vx / y1.uR.rho);
            double vyS = 0.5 * (x.uL.vy / x.uL.rho + x.uR.vy / x.uR.rho);
            double vyN = 0.5 * (x1.uL.vy / x1.uL.rho + x1.uR.vy / x1.uR.rho);

            double ByW = Cn.Byf[j + nghost][i + nghost - 1];
            double ByE = Cn.Byf[j + nghost][i + nghost];
            double BxS = Cn.Bxf[j + nghost - 1][i + nghost];
            double BxN = Cn.Bxf[j + nghost][i + nghost];

            Ez[j][i] = - (aW*vxW*ByW + aE*vxE*ByE) + (aS*vyS*BxS + aN*vyN*BxN) + (dE*ByE - dW*ByW) - (dN*BxN - dS*BxS);
        }
    }

    return Ez;

}

#include "WrittingUtils.hpp"

void ApplyConstrainedTransport(ConservativeVariables& Cn1, const ConservativeVariables& Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, Slope sl, Riemann rs, CTMethod ct)
{
    CTFunction ChosenCT = getCT(ct);
    std::vector<std::vector<double>> Ez = ChosenCT(Cn, Dx, Dy, Dt, nghost, rec, sl, rs);
    
    for(int j = nghost; j < Cn1.ny - nghost; ++j)
    {
        for(int i = nghost; i <= Cn1.nx - nghost; ++i)
        {
            Cn1.Bxf[j][i] = Cn.Bxf[j][i] - (Dt / Dy) * (Ez[j + 1 - nghost][i - nghost] - Ez[j - nghost][i - nghost]);
        }
    }

    for(int j = nghost; j <= Cn1.ny - nghost; ++j)
    {
        for(int i = nghost; i < Cn1.nx - nghost; ++i)
        {
            Cn1.Byf[j][i] = Cn.Byf[j][i] + (Dt / Dx) * (Ez[j - nghost][i + 1 - nghost] - Ez[j - nghost][i - nghost]);
        }
    }
}
