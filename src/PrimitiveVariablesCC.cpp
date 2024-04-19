#include "ConservativeVariablesCC.hpp"
#include "PrimitiveVariablesCC.hpp"


PrimitiveVariablesCC::PrimitiveVariablesCC(int nx_, int ny_) : nx(nx_), ny(ny_)
{
    rho.resize(nx, std::vector<double>(ny));
    vx.resize(nx, std::vector<double>(ny));
    vy.resize(nx, std::vector<double>(ny));
    vz.resize(nx, std::vector<double>(ny));
    Bx.resize(nx, std::vector<double>(ny));
    By.resize(nx, std::vector<double>(ny));
    Bz.resize(nx, std::vector<double>(ny));
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

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            rho[i][j] = C_cc.rho[i][j];
            vx[i][j] = C_cc.rhovx[i][j] / C_cc.rho[i][j];
            vy[i][j] = C_cc.rhovy[i][j] / C_cc.rho[i][j];
            vz[i][j] = C_cc.rhovz[i][j] / C_cc.rho[i][j];
            Bx[i][j] = C_cc.Bx[i][j];
            By[i][j] = C_cc.By[i][j];
            Bz[i][j] = C_cc.Bz[i][j];
        }
    }
}

PrimitiveVariablesCC::~PrimitiveVariablesCC() = default;

void PrimitiveVariablesCC::set(const ReconstructedValues& rv, int i, int j){
    rho[i][j] = rv.rho;
    vx[i][j] = rv.vx;
    vy[i][j] = rv.vy;
    vz[i][j] = rv.vz;
    Bx[i][j] = rv.Bx;
    By[i][j] = rv.By;
    Bz[i][j] = rv.Bz;
}

void PrimitiveVariablesCC::init(const std::vector<std::vector<double>>& rho_, const std::vector<std::vector<double>>& vx_, const std::vector<std::vector<double>>& vy_, const std::vector<std::vector<double>>& vz_, const std::vector<std::vector<double>>& Bx_, const std::vector<std::vector<double>>& By_, const std::vector<std::vector<double>>& Bz_, const double P_){
    rho = rho_;
    vx = vx_;
    vy = vy_;
    vz = vz_;
    Bx = Bx_;
    By = By_;
    Bz = Bz_;
}

ReconstructedValues PrimitiveVariablesCC::operator()(int i, int j) const {
    if (i < 0 || i >= nx || j < 0 || j >= ny)
        throw("Index out of range");

    return ReconstructedValues{rho[i][j], vx[i][j], vy[i][j], vz[i][j], Bx[i][j], By[i][j], Bz[i][j]};
}
