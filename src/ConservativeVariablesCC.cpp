#include "PrimitiveVariablesCC.hpp"
#include "ConservativeVariablesCC.hpp"

ConservativeVariablesCC::ConservativeVariablesCC(int nx, int ny) : nx(nx), ny(ny)
{
    rho.resize(nx, std::vector<double>(ny));
    rhovx.resize(nx, std::vector<double>(ny));
    rhovy.resize(nx, std::vector<double>(ny));
    rhovz.resize(nx, std::vector<double>(ny));
    Bx.resize(nx, std::vector<double>(ny));
    By.resize(nx, std::vector<double>(ny));
    Bz.resize(nx, std::vector<double>(ny));
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

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            rho[j][i] = P_cc.rho[j][i];
            rhovx[j][i] = P_cc.rho[j][i] * P_cc.vx[j][i];
            rhovy[j][i] = P_cc.rho[j][i] * P_cc.vy[j][i];
            rhovz[j][i] = P_cc.rho[j][i] * P_cc.vz[j][i];
            Bx[j][i] = P_cc.Bx[j][i];
            By[j][i] = P_cc.By[j][i];
            Bz[j][i] = P_cc.Bz[j][i];
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
}

ReconstructedValues ConservativeVariablesCC::operator()(int i, int j) const {
    if (i < 0 || i >= nx || j < 0 || j >= ny)
        throw("Index out of range");

    return ReconstructedValues{rho[j][i], rhovx[j][i], rhovy[j][i], rhovz[j][i], Bx[j][i], By[j][i], Bz[j][i]};
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
        }
    }
    return result;
}

ConservativeVariablesCC ConservativeVariablesCC::operator-(const ConservativeVariablesCC& other) const {
    ConservativeVariablesCC result(this->nx, this->ny);

    if (this->nx != other.nx || this->ny != other.ny) {
        throw("Dimension mismatch");
    }

    for (int i = 0; i < this->nx; ++i) {
        for (int j = 0; j < this->ny; ++j) {
            result.rho[i][j] = this->rho[i][j] - other.rho[i][j];
            result.rhovx[i][j] = this->rhovx[i][j] - other.rhovx[i][j];
            result.rhovy[i][j] = this->rhovy[i][j] - other.rhovy[i][j];
            result.rhovz[i][j] = this->rhovz[i][j] - other.rhovz[i][j];
            result.Bx[i][j] = this->Bx[i][j] - other.Bx[i][j];
            result.By[i][j] = this->By[i][j] - other.By[i][j];
            result.Bz[i][j] = this->Bz[i][j] - other.Bz[i][j];
        }
    }
    return result;
}

ConservativeVariablesCC ConservativeVariablesCC::operator+(const ConservativeVariablesCC& other) const {
    ConservativeVariablesCC result(this->nx, this->ny);

    if (this->nx != other.nx || this->ny != other.ny) {
        throw("Dimension mismatch");
    }

    for (int i = 0; i < this->nx; ++i) {
        for (int j = 0; j < this->ny; ++j) {
            result.rho[i][j] = this->rho[i][j] + other.rho[i][j];
            result.rhovx[i][j] = this->rhovx[i][j] + other.rhovx[i][j];
            result.rhovy[i][j] = this->rhovy[i][j] + other.rhovy[i][j];
            result.rhovz[i][j] = this->rhovz[i][j] + other.rhovz[i][j];
            result.Bx[i][j] = this->Bx[i][j] + other.Bx[i][j];
            result.By[i][j] = this->By[i][j] + other.By[i][j];
            result.Bz[i][j] = this->Bz[i][j] + other.Bz[i][j];
        }
    }
    return result;
}

ConservativeVariablesCC& ConservativeVariablesCC::operator=(const ConservativeVariablesCC& other) {
    if (this->nx != other.nx || this->ny != other.ny) {
        throw("Dimension mismatch");
    }

    for (int i = 0; i < this->nx; ++i) {
        for (int j = 0; j < this->ny; ++j) {
            this->rho[i][j] = other.rho[i][j];
            this->rhovx[i][j] = other.rhovx[i][j];
            this->rhovy[i][j] = other.rhovy[i][j];
            this->rhovz[i][j] = other.rhovz[i][j];
            this->Bx[i][j] = other.Bx[i][j];
            this->By[i][j] = other.By[i][j];
            this->Bz[i][j] = other.Bz[i][j];
        }
    }
    return *this;
}
