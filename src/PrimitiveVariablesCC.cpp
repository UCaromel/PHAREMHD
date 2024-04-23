#include "ConservativeVariablesCC.hpp"
#include "PrimitiveVariablesCC.hpp"
#include "EquationOfState.hpp"


PrimitiveVariablesCC::PrimitiveVariablesCC(int nx_, int ny_) : nx(nx_), ny(ny_)
{
    rho.resize(nx, std::vector<double>(ny));
    vx.resize(nx, std::vector<double>(ny));
    vy.resize(nx, std::vector<double>(ny));
    vz.resize(nx, std::vector<double>(ny));
    Bx.resize(nx, std::vector<double>(ny));
    By.resize(nx, std::vector<double>(ny));
    Bz.resize(nx, std::vector<double>(ny));
    P.resize(nx, std::vector<double>(ny));
}

PrimitiveVariablesCC::PrimitiveVariablesCC(const ConservativeVariablesCC& C_cc) : nx(C_cc.nx), ny(C_cc.ny)
{
    rho.resize(nx, std::vector<double>(ny));
    vx.resize(nx, std::vector<double>(ny));
    vy.resize(nx, std::vector<double>(ny));
    vz.resize(nx, std::vector<double>(ny));
    Bx.resize(nx, std::vector<double>(ny));
    By.resize(nx, std::vector<double>(ny));
    Bz.resize(nx, std::vector<double>(ny));
    P.resize(nx, std::vector<double>(ny));

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
        }
    }
}

PrimitiveVariablesCC::~PrimitiveVariablesCC() = default;

void PrimitiveVariablesCC::set(const ReconstructedValues& rv, int i, int j){
    rho[j][i] = rv.rho;
    vx[j][i] = rv.vx;
    vy[j][i] = rv.vy;
    vz[j][i] = rv.vz;
    Bx[j][i] = rv.Bx;
    By[j][i] = rv.By;
    Bz[j][i] = rv.Bz;
    P[j][i] = rv.P;
}

void PrimitiveVariablesCC::init(const std::vector<std::vector<double>>& rho_, const std::vector<std::vector<double>>& vx_, const std::vector<std::vector<double>>& vy_, const std::vector<std::vector<double>>& vz_, const std::vector<std::vector<double>>& Bx_, const std::vector<std::vector<double>>& By_, const std::vector<std::vector<double>>& Bz_, const std::vector<std::vector<double>>& P_){
    rho = rho_;
    vx = vx_;
    vy = vy_;
    vz = vz_;
    Bx = Bx_;
    By = By_;
    Bz = Bz_;
    P = P_;
}

ReconstructedValues PrimitiveVariablesCC::operator()(int i, int j) const {
    if (i < 0 || i >= nx || j < 0 || j >= ny)
        throw("Index out of range");

    return ReconstructedValues{rho[j][i], vx[j][i], vy[j][i], vz[j][i], Bx[j][i], By[j][i], Bz[j][i], P[j][i]};
}

PrimitiveVariablesCC& PrimitiveVariablesCC::operator=(const PrimitiveVariablesCC& other) {
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
        }
    }
    return *this;
}