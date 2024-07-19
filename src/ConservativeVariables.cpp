#include "PrimitiveVariables.hpp"
#include "ConservativeVariables.hpp"
#include "EquationOfState.hpp"

ConservativeVariables::ConservativeVariables(int nx, int ny) : nx(nx), ny(ny)
{
    rho.resize(ny, std::vector<double>(nx));
    rhovx.resize(ny, std::vector<double>(nx));
    rhovy.resize(ny, std::vector<double>(nx));
    rhovz.resize(ny, std::vector<double>(nx));
    Bx.resize(ny, std::vector<double>(nx));
    By.resize(ny, std::vector<double>(nx));
    Bz.resize(ny, std::vector<double>(nx));
    Etot.resize(ny, std::vector<double>(nx));

    Bxf.resize(ny, std::vector<double>(nx + 1));
    Byf.resize(ny + 1, std::vector<double>(nx));

    // J needs one extra ghost cell
    int nxJ = nx + 2;
    int nyJ = ny + 2;
    Jx.resize(nyJ + 1, std::vector<double>(nxJ));
    Jy.resize(nyJ, std::vector<double>(nxJ + 1));
    Jz.resize(nyJ + 1, std::vector<double>(nxJ + 1));
}

ConservativeVariables::ConservativeVariables(const PrimitiveVariables& P_cc) : nx(P_cc.nx), ny(P_cc.ny)
{
    rho.resize(ny, std::vector<double>(nx));
    rhovx.resize(ny, std::vector<double>(nx));
    rhovy.resize(ny, std::vector<double>(nx));
    rhovz.resize(ny, std::vector<double>(nx));
    Bx.resize(ny, std::vector<double>(nx));
    By.resize(ny, std::vector<double>(nx));
    Bz.resize(ny, std::vector<double>(nx));
    Etot.resize(ny, std::vector<double>(nx));

    Bxf.resize(ny, std::vector<double>(nx + 1));
    Byf.resize(ny + 1, std::vector<double>(nx));

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            rho[j][i] = P_cc.rho[j][i];
            rhovx[j][i] = P_cc.rho[j][i] * P_cc.vx[j][i];
            rhovy[j][i] = P_cc.rho[j][i] * P_cc.vy[j][i];
            rhovz[j][i] = P_cc.rho[j][i] * P_cc.vz[j][i];
            Bx[j][i] = P_cc.Bx[j][i];
            By[j][i] = P_cc.By[j][i];
            Bz[j][i] = P_cc.Bz[j][i];
            Etot[j][i] = EosEtot(P_cc(i,j));

            Bxf[j][i] = P_cc.Bxf[j][i];
            Byf[j][i] = P_cc.Byf[j][i];
        }
        Bxf[j][nx] = P_cc.Bxf[j][nx];
    }
    for (int i = 0; i < nx; i++) {
        Byf[ny][i] = P_cc.Byf[ny][i];
    }

    // J needs one extra ghost cell
    int nxJ = nx + 2;
    int nyJ = ny + 2;
    Jx.resize(nyJ + 1, std::vector<double>(nxJ));
    Jy.resize(nyJ, std::vector<double>(nxJ + 1));
    Jz.resize(nyJ + 1, std::vector<double>(nxJ + 1));

    for (int j = 0; j < nyJ; j++) {
        for (int i = 0; i < nxJ; i++) {
            Jx[j][i] = P_cc.Jx[j][i];
            Jy[j][i] = P_cc.Jy[j][i];
            Jz[j][i] = P_cc.Jz[j][i];
        }
        Jy[j][nxJ] = P_cc.Jx[j][nxJ];
        Jz[j][nxJ] = P_cc.Jx[j][nxJ];
    }
    for (int i = 0; i < nxJ; i++) {
        Jx[nyJ][i] = P_cc.Jx[nyJ][i];
        Jz[nyJ][i] = P_cc.Jx[nyJ][i];
    }
    Jz[nyJ][nxJ] = P_cc.Jx[nyJ][nxJ];
}

ConservativeVariables::~ConservativeVariables() = default;

void ConservativeVariables::ReconstructCenteredB(int nghost) {
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

void ConservativeVariables::set(const ReconstructedValues& rv, int i, int j){
    rho[j][i] = rv.rho;
    rhovx[j][i] = rv.vx;
    rhovy[j][i] = rv.vy;
    rhovz[j][i] = rv.vz;
    Bx[j][i] = rv.Bx;
    By[j][i] = rv.By;
    Bz[j][i] = rv.Bz;
    Etot[j][i] = rv.P;
}

ReconstructedValues ConservativeVariables::operator()(int i, int j) const {
    if (i < 0 || i >= nx || j < 0 || j >= ny)
        throw("Index out of range");

    return ReconstructedValues{rho[j][i], rhovx[j][i], rhovy[j][i], rhovz[j][i], Bx[j][i], By[j][i], Bz[j][i], Etot[j][i]};
}


ConservativeVariables ConservativeVariables::operator*(double scalar) const {
    ConservativeVariables result(this->nx, this->ny);

    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            result.rho[j][i] = this->rho[j][i] * scalar;
            result.rhovx[j][i] = this->rhovx[j][i] * scalar;
            result.rhovy[j][i] = this->rhovy[j][i] * scalar;
            result.rhovz[j][i] = this->rhovz[j][i] * scalar;
            result.Bx[j][i] = this->Bx[j][i] * scalar;
            result.By[j][i] = this->By[j][i] * scalar;
            result.Bz[j][i] = this->Bz[j][i] * scalar;
            result.Etot[j][i] = this->Etot[j][i] * scalar;

            result.Bxf[j][i] = this->Bxf[j][i] * scalar;
            result.Byf[j][i] = this->Byf[j][i] * scalar;
        }
        result.Bxf[j][nx] = this->Bxf[j][nx] * scalar;
    }
    for (int i = 0; i < nx; i++) {
        result.Byf[ny][i] = this->Byf[ny][i] * scalar;
    }
    return result;
}

ConservativeVariables ConservativeVariables::operator-(const ConservativeVariables& other) const {
    ConservativeVariables result(this->nx, this->ny);

    if (this->nx != other.nx || this->ny != other.ny) {
        throw("Dimension mismatch");
    }

    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            result.rho[j][i] = this->rho[j][i] - other.rho[j][i];
            result.rhovx[j][i] = this->rhovx[j][i] - other.rhovx[j][i];
            result.rhovy[j][i] = this->rhovy[j][i] - other.rhovy[j][i];
            result.rhovz[j][i] = this->rhovz[j][i] - other.rhovz[j][i];
            result.Bx[j][i] = this->Bx[j][i] - other.Bx[j][i];
            result.By[j][i] = this->By[j][i] - other.By[j][i];
            result.Bz[j][i] = this->Bz[j][i] - other.Bz[j][i];
            result.Etot[j][i] = this->Etot[j][i] - other.Etot[j][i];

            result.Bxf[j][i] = this->Bxf[j][i] - other.Bxf[j][i];
            result.Byf[j][i] = this->Byf[j][i] - other.Byf[j][i];
        }
        result.Bxf[j][nx] = this->Bxf[j][nx] - other.Bxf[j][nx];
    }
    for (int i = 0; i < nx; i++) {
        result.Byf[ny][i] = this->Byf[ny][i] - other.Byf[ny][i];
    }
    return result;
}

ConservativeVariables ConservativeVariables::operator+(const ConservativeVariables& other) const {
    ConservativeVariables result(this->nx, this->ny);

    if (this->nx != other.nx || this->ny != other.ny) {
        throw("Dimension mismatch");
    }

    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            result.rho[j][i] = this->rho[j][i] + other.rho[j][i];
            result.rhovx[j][i] = this->rhovx[j][i] + other.rhovx[j][i];
            result.rhovy[j][i] = this->rhovy[j][i] + other.rhovy[j][i];
            result.rhovz[j][i] = this->rhovz[j][i] + other.rhovz[j][i];
            result.Bx[j][i] = this->Bx[j][i] + other.Bx[j][i];
            result.By[j][i] = this->By[j][i] + other.By[j][i];
            result.Bz[j][i] = this->Bz[j][i] + other.Bz[j][i];
            result.Etot[j][i] = this->Etot[j][i] + other.Etot[j][i];

            result.Bxf[j][i] = this->Bxf[j][i] + other.Bxf[j][i];
            result.Byf[j][i] = this->Byf[j][i] + other.Byf[j][i];
        }
        result.Bxf[j][nx] = this->Bxf[j][nx] + other.Bxf[j][nx];
    }
    for (int i = 0; i < nx; i++) {
        result.Byf[ny][i] = this->Byf[ny][i] + other.Byf[ny][i];
    }
    return result;
}

ConservativeVariables& ConservativeVariables::operator=(const ConservativeVariables& other) {
    if (this->nx != other.nx || this->ny != other.ny) {
        throw("Dimension mismatch");
    }

    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            this->rho[j][i] = other.rho[j][i];
            this->rhovx[j][i] = other.rhovx[j][i];
            this->rhovy[j][i] = other.rhovy[j][i];
            this->rhovz[j][i] = other.rhovz[j][i];
            this->Bx[j][i] = other.Bx[j][i];
            this->By[j][i] = other.By[j][i];
            this->Bz[j][i] = other.Bz[j][i];
            this->Etot[j][i] = other.Etot[j][i];

            this->Bxf[j][i] = other.Bxf[j][i];
            this->Byf[j][i] = other.Byf[j][i];
        }
        this->Bxf[j][nx] = other.Bxf[j][nx];
    }
    for (int i = 0; i < nx; i++) {
        this->Byf[ny][i] = other.Byf[ny][i];
    }

    return *this;
}
