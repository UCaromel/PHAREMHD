#include "PrimitiveVariablesCC.hpp"
#include "ConservativeVariablesCC.hpp"
#include "EquationOfState.hpp"

ConservativeVariablesCC::ConservativeVariablesCC(int nx, int ny) : nx(nx), ny(ny)
{
    rho.resize(nx, std::vector<double>(ny));
    rhovx.resize(nx, std::vector<double>(ny));
    rhovy.resize(nx, std::vector<double>(ny));
    rhovz.resize(nx, std::vector<double>(ny));
    Bx.resize(nx, std::vector<double>(ny));
    By.resize(nx, std::vector<double>(ny));
    Bz.resize(nx, std::vector<double>(ny));
    Etot.resize(nx, std::vector<double>(ny));
}

ConservativeVariablesCC::ConservativeVariablesCC(const PrimitiveVariablesCC& P_cc) : nx(P_cc.nx), ny(P_cc.ny)
{
    rho.resize(nx, std::vector<double>(ny));
    rhovx.resize(nx, std::vector<double>(ny));
    rhovy.resize(nx, std::vector<double>(ny));
    rhovz.resize(nx, std::vector<double>(ny));
    Bx.resize(nx, std::vector<double>(ny));
    By.resize(nx, std::vector<double>(ny));
    Bz.resize(nx, std::vector<double>(ny));
    Etot.resize(nx, std::vector<double>(ny));

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
        }
    }
}

ConservativeVariablesCC::~ConservativeVariablesCC() = default;

void ConservativeVariablesCC::set(const ReconstructedValues& rv, int i, int j){
    rho[j][i] = rv.rho;
    rhovx[j][i] = rv.vx;
    rhovy[j][i] = rv.vy;
    rhovz[j][i] = rv.vz;
    Bx[j][i] = rv.Bx;
    By[j][i] = rv.By;
    Bz[j][i] = rv.Bz;
    Etot[j][i] = rv.P;
}

ReconstructedValues ConservativeVariablesCC::operator()(int i, int j) const {
    if (i < 0 || i >= nx || j < 0 || j >= ny)
        throw("Index out of range");

    return ReconstructedValues{rho[j][i], rhovx[j][i], rhovy[j][i], rhovz[j][i], Bx[j][i], By[j][i], Bz[j][i], Etot[j][i]};
}

ConservativeVariablesCC ConservativeVariablesCC::operator*(double scalar) const {
    ConservativeVariablesCC result(this->nx, this->ny);

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
        }
    }
    return result;
}

ConservativeVariablesCC ConservativeVariablesCC::operator-(const ConservativeVariablesCC& other) const {
    ConservativeVariablesCC result(this->nx, this->ny);

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
        }
    }
    return result;
}

ConservativeVariablesCC ConservativeVariablesCC::operator+(const ConservativeVariablesCC& other) const {
    ConservativeVariablesCC result(this->nx, this->ny);

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
        }
    }
    return result;
}

ConservativeVariablesCC& ConservativeVariablesCC::operator=(const ConservativeVariablesCC& other) {
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
        }
    }
    return *this;
}
