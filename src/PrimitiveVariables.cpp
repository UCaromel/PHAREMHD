#include "ConservativeVariables.hpp"
#include "PrimitiveVariables.hpp"
#include "EquationOfState.hpp"


PrimitiveVariables::PrimitiveVariables(int nx_, int ny_) : nx(nx_), ny(ny_)
{
    rho.resize(ny, std::vector<double>(nx));
    vx.resize(ny, std::vector<double>(nx));
    vy.resize(ny, std::vector<double>(nx));
    vz.resize(ny, std::vector<double>(nx));
    Bx.resize(ny, std::vector<double>(nx));
    By.resize(ny, std::vector<double>(nx));
    Bz.resize(ny, std::vector<double>(nx));
    P.resize(ny, std::vector<double>(nx));

    Bxf.resize(ny, std::vector<double>(nx + 1));
    Byf.resize(ny + 1, std::vector<double>(nx));

    // J needs one extra ghost cell
    int nxJ = nx + 2;
    int nyJ = ny + 2;
    Jx.resize(nyJ + 1, std::vector<double>(nxJ));
    Jy.resize(nyJ, std::vector<double>(nxJ + 1));
    Jz.resize(nyJ + 1, std::vector<double>(nxJ + 1));
}

PrimitiveVariables::PrimitiveVariables(const ConservativeVariables& C_cc) : nx(C_cc.nx), ny(C_cc.ny)
{
    rho.resize(ny, std::vector<double>(nx));
    vx.resize(ny, std::vector<double>(nx));
    vy.resize(ny, std::vector<double>(nx));
    vz.resize(ny, std::vector<double>(nx));
    Bx.resize(ny, std::vector<double>(nx));
    By.resize(ny, std::vector<double>(nx));
    Bz.resize(ny, std::vector<double>(nx));
    P.resize(ny, std::vector<double>(nx));

    Bxf.resize(ny, std::vector<double>(nx + 1));
    Byf.resize(ny + 1, std::vector<double>(nx));

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            rho[j][i] = C_cc.rho[j][i];
            vx[j][i] = C_cc.rhovx[j][i] / C_cc.rho[j][i];
            vy[j][i] = C_cc.rhovy[j][i] / C_cc.rho[j][i];
            vz[j][i] = C_cc.rhovz[j][i] / C_cc.rho[j][i];
            Bx[j][i] = C_cc.Bx[j][i];
            By[j][i] = C_cc.By[j][i];
            Bz[j][i] = C_cc.Bz[j][i];
            P[j][i] = EosP(C_cc(i,j));

            Bxf[j][i] = C_cc.Bxf[j][i];
            Byf[j][i] = C_cc.Byf[j][i];
        }
        Bxf[j][nx] = C_cc.Bxf[j][nx];
    }
    for (int i = 0; i < nx; i++) {
        Byf[ny][i] = C_cc.Byf[ny][i];
    }


    // J needs one extra ghost cell
    int nxJ = nx + 2;
    int nyJ = ny + 2;
    Jx.resize(nyJ + 1, std::vector<double>(nxJ));
    Jy.resize(nyJ, std::vector<double>(nxJ + 1));
    Jz.resize(nyJ + 1, std::vector<double>(nxJ + 1));

    for (int j = 0; j < nyJ; j++) {
        for (int i = 0; i < nxJ; i++) {
            Jx[j][i] = C_cc.Jx[j][i];
            Jy[j][i] = C_cc.Jy[j][i];
            Jz[j][i] = C_cc.Jz[j][i];
        }
        Jy[j][nxJ] = C_cc.Jx[j][nxJ];
        Jz[j][nxJ] = C_cc.Jx[j][nxJ];
    }
    for (int i = 0; i < nxJ; i++) {
        Jx[nyJ][i] = C_cc.Jx[nyJ][i];
        Jz[nyJ][i] = C_cc.Jx[nyJ][i];
    }
    Jz[nyJ][nxJ] = C_cc.Jx[nyJ][nxJ];

}

PrimitiveVariables::~PrimitiveVariables() = default;

void PrimitiveVariables::ReconstructCenteredB(int nghost) {
    for (int i = 0; i < nx - 2*nghost; i++) {
        for (int j = 0; j < ny - 2*nghost; j++) {
            if (nghost == 0) {
                Bx[j][i] = 0.5*(Bxf[j][i] + Bxf[j][i + 1]);
                By[j][i] = 0.5*(Byf[j][i] + Byf[j + 1][i]);
            } else {
                Bx[j + nghost][i + nghost] = 0.5*(Bxf[j + nghost][i + nghost] + Bxf[j + nghost][i + nghost + 1]);
                By[j + nghost][i + nghost] = 0.5*(Byf[j + nghost][i + nghost] + Byf[j + nghost + 1][i + nghost]);
            }
        }    
    }
}

void PrimitiveVariables::init(const std::vector<std::vector<double>>& rho_, const std::vector<std::vector<double>>& vx_, const std::vector<std::vector<double>>& vy_, const std::vector<std::vector<double>>& vz_, const std::vector<std::vector<double>>& Bxf_, const std::vector<std::vector<double>>& Byf_, const std::vector<std::vector<double>>& Bz_, const std::vector<std::vector<double>>& P_){
    rho = rho_;
    vx = vx_;
    vy = vy_;
    vz = vz_;
    Bxf = Bxf_;
    Byf = Byf_;
    Bz = Bz_;
    P = P_;

    this -> ReconstructCenteredB(0);
}

void PrimitiveVariables::set(const ReconstructedValues& rv, int i, int j){
    rho[j][i] = rv.rho;
    vx[j][i] = rv.vx;
    vy[j][i] = rv.vy;
    vz[j][i] = rv.vz;
    Bx[j][i] = rv.Bx;
    By[j][i] = rv.By;
    Bz[j][i] = rv.Bz;
    P[j][i] = rv.P;
}

ReconstructedValues PrimitiveVariables::operator()(int i, int j) const {
    if (i < 0 || i >= nx || j < 0 || j >= ny)
        throw("Index out of range");

    return ReconstructedValues{rho[j][i], vx[j][i], vy[j][i], vz[j][i], Bx[j][i], By[j][i], Bz[j][i], P[j][i]};
}

PrimitiveVariables& PrimitiveVariables::operator=(const PrimitiveVariables& other) {
    if (this->nx != other.nx || this->ny != other.ny) {
        throw("Dimension mismatch");
    }

    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            this->rho[j][i] = other.rho[j][i];
            this->vx[j][i] = other.vx[j][i];
            this->vy[j][i] = other.vy[j][i];
            this->vz[j][i] = other.vz[j][i];
            this->Bx[j][i] = other.Bx[j][i];
            this->By[j][i] = other.By[j][i];
            this->Bz[j][i] = other.Bz[j][i];
            this->P[j][i] = other.P[j][i];

            this->Bxf[j][i] = other.Bxf[j][i];
            this->Byf[j][i] = other.Byf[j][i];
        }
        Bxf[j][nx] = other.Bxf[j][nx];
    }
    for (int i = 0; i < nx; i++) {
        Byf[ny][i] = other.Byf[ny][i];
    }

    // J needs one extra ghost cell
    int nxJ = nx + 2;
    int nyJ = ny + 2;

    for (int j = 0; j < nyJ; j++) {
        for (int i = 0; i < nxJ; i++) {
            this->Jx[j][i] = other.Jx[j][i];
            this->Jy[j][i] = other.Jy[j][i];
            this->Jz[j][i] = other.Jz[j][i];
        }
        this->Jy[j][nxJ] = other.Jx[j][nxJ];
        this->Jz[j][nxJ] = other.Jx[j][nxJ];
    }
    for (int i = 0; i < nxJ; i++) {
        this->Jx[nyJ][i] = other.Jx[nyJ][i];
        this->Jz[nyJ][i] = other.Jx[nyJ][i];
    }
    this->Jz[nyJ][nxJ] = other.Jx[nyJ][nxJ];

    return *this;
}