#include "ConstainedTransport.hpp"

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> ConstrainedTransportAverage(const ConservativeVariablesCC &Cn /* Assuming ghost cells are added */, double Dx, double Dy, double Dt, int nghost, Reconstruction rec)
{
    // Edge-centered
    std::vector<std::vector<double>> vx(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
    std::vector<std::vector<double>> vy(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
    std::vector<std::vector<double>> Bx(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
    std::vector<std::vector<double>> By(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
    std::vector<std::vector<double>> Ez(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));

    // Face-centered
    std::vector<std::vector<double>> bx(Cn.ny - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
    std::vector<std::vector<double>> by(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx - 2 * nghost));

    std::vector<std::vector<double>> BX(Cn.ny - 2 * nghost, std::vector<double>(Cn.nx - 2 * nghost));
    std::vector<std::vector<double>> BY(Cn.ny - 2 * nghost, std::vector<double>(Cn.nx - 2 * nghost));

    for (int j = nghost; j <= Cn.ny - nghost; ++j)
    {
        for (int i = nghost; i <= Cn.nx - nghost; ++i)
        {
            vx[j - nghost][i - nghost] = 0.25 * ((Cn.rhovx[j][i] / Cn.rho[j][i]) + (Cn.rhovx[j][i - 1] / Cn.rho[j][i - 1]) + (Cn.rhovx[j - 1][i] / Cn.rho[j - 1][i]) + (Cn.rhovx[j - 1][i - 1] / Cn.rho[j - 1][i - 1]));
            vy[j - nghost][i - nghost] = 0.25 * ((Cn.rhovy[j][i] / Cn.rho[j][i]) + (Cn.rhovy[j][i - 1] / Cn.rho[j][i - 1]) + (Cn.rhovy[j - 1][i] / Cn.rho[j - 1][i]) + (Cn.rhovy[j - 1][i - 1] / Cn.rho[j - 1][i - 1]));
            Bx[j - nghost][i - nghost] = 0.25 * (Cn.Bx[j][i] + Cn.Bx[j][i - 1] + Cn.Bx[j - 1][i] + Cn.Bx[j - 1][i - 1]);
            By[j - nghost][i - nghost] = 0.25 * (Cn.By[j][i] + Cn.By[j][i - 1] + Cn.By[j - 1][i] + Cn.By[j - 1][i - 1]);
            Ez[j - nghost][i - nghost] = vy[j - nghost][i - nghost] * Bx[j - nghost][i - nghost] - vx[j - nghost][i - nghost] * By[j - nghost][i - nghost];
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

    return {BX, BY};
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> UCTHLL(const ConservativeVariablesCC &Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec){
    PrimitiveVariablesCC Pn(Cn);
    std::vector<std::vector<Interface>> InterfacesX(Pn.ny + 2 - 2 * nghost, std::vector<Interface>(Pn.nx + 1 - 2 * nghost));
    std::vector<std::vector<Interface>> InterfacesY(Pn.ny + 1 - 2 * nghost, std::vector<Interface>(Pn.nx + 2 - 2 * nghost));

    std::vector<std::vector<double>> bx(Cn.ny - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
    std::vector<std::vector<double>> by(Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx - 2 * nghost));

    std::vector<std::vector<double>> BX(Cn.ny - 2 * nghost, std::vector<double>(Cn.nx - 2 * nghost));
    std::vector<std::vector<double>> BY(Cn.ny - 2 * nghost, std::vector<double>(Cn.nx - 2 * nghost));

    std::vector<std::vector<double>> Ez(Pn.ny + 1 - 2 * nghost, std::vector<double>(Pn.nx + 1 - 2 * nghost));

    for (int j = nghost - 1; j < Cn.ny - nghost + 1; ++j)
    {
        for (int i = nghost; i <= Cn.nx - nghost; ++i)
        {
            InterfacesX[j - nghost + 1][i - nghost] = Interface(Pn, i, j, rec, nghost, Dir::X);
        }
    }

    for (int j = nghost; j <= Cn.ny - nghost; ++j)
    {
        for (int i = nghost - 1; i < Cn.nx - nghost + 1; ++i)
        {
            InterfacesY[j - nghost][i - nghost + 1] = Interface(Pn, i, j, rec, nghost, Dir::Y);
        }
    }

    for (int j = 0; j <= Cn.ny - 2*nghost; ++j){
        for (int i = 0; i <= Cn.nx - 2*nghost; ++i){
            Interface x = InterfacesX[j][i];
            Interface x1 = InterfacesX[j + 1][i];
            Interface y = InterfacesY[j][i];
            Interface y1 = InterfacesY[j][i + 1];

            // Compute solver's coeficients
            double alphaxR = std::max(0.0, x.uR.vx + x.cfastxR);
            double alphaxL = - std::min(0.0, x.uL.vx - x.cfastxL);

            double alphax1R = std::max(0.0, x1.uR.vx + x1.cfastxR);
            double alphax1L = - std::min(0.0, x1.uL.vx - x1.cfastxL);

            double alphayR = std::max(0.0, y.uR.vy + y.cfastyR);
            double alphayL = - std::min(0.0, y.uL.vy - y.cfastyL);

            double alphay1R = std::max(0.0, y1.uR.vy + y1.cfastyR);
            double alphay1L = - std::min(0.0, y1.uL.vy - y1.cfastyL);


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

            double vxW = 0.5 * (y.uL.vx + y.uR.vx);
            double vxE = 0.5 * (y1.uL.vx + y1.uR.vx);
            double vyS = 0.5 * (x.uL.vy + x.uR.vy);
            double vyN = 0.5 * (x1.uL.vy + x1.uR.vy);

            double ByW = 0.5 * (y.uL.By + y.uR.By);
            double ByE = 0.5 * (y1.uL.By + y1.uR.By);
            double BxS = 0.5 * (x.uL.Bx + x.uR.Bx);
            double BxN = 0.5 * (x1.uL.Bx + x1.uR.Bx);

            Ez[j][i] = - (aW*vxW*ByW + aE*vxE*ByE) + (aS*vyS*BxS + aN*vyN*BxN) + (dE*ByE - dW*ByW) - (dN*BxN - dS*BxS);
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

    return {BX, BY};

}

void ApplyConstrainedTransport(ConservativeVariablesCC& Cn1, const ConservativeVariablesCC& Cn, double Dx, double Dy, double Dt, int nghost, Reconstruction rec, CTMethod ct)
{
    CTFunction ChosenCT = getCT(ct);
    auto [BX, BY] = ChosenCT(Cn, Dx, Dy, Dt, nghost, rec);
    
    for(int j = nghost; j < Cn1.ny - nghost; ++j)
    {
        for(int i = nghost; i < Cn1.nx - nghost; ++i)
        {
            Cn1.Bx[j][i] = BX[j - nghost][i - nghost];
            Cn1.By[j][i] = BY[j - nghost][i - nghost];
        }
    }
}
