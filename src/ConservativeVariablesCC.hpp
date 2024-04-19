#ifndef CONSERVATIVE_VARIABLES_CC_HPP_
#define CONSERVATIVE_VARIABLES_CC_HPP_

#include <vector>
#include <string>

#include "PrimitiveVariablesCC.hpp"

class PrimitiveVariablesCC;

class ConservativeVariablesCC {
public:
    int nx;
    int ny;
    std::vector<std::vector<double>> rho;
    std::vector<std::vector<double>> rhovx;
    std::vector<std::vector<double>> rhovy;
    std::vector<std::vector<double>> rhovz;
    std::vector<std::vector<double>> Bx;
    std::vector<std::vector<double>> By;
    std::vector<std::vector<double>> Bz;
    double P;

        // Constructor that takes grid dimensions as arguments
    ConservativeVariablesCC(int nx, int ny) : nx(nx), ny(ny)
    {
        rho.resize(nx, std::vector<double>(ny));
        rhovx.resize(nx, std::vector<double>(ny));
        rhovy.resize(nx, std::vector<double>(ny));
        rhovz.resize(nx, std::vector<double>(ny));
        Bx.resize(nx, std::vector<double>(ny));
        By.resize(nx, std::vector<double>(ny));
        Bz.resize(nx, std::vector<double>(ny));
    }

    ConservativeVariablesCC(const PrimitiveVariablesCC& P_cc) : nx(P_cc.nx), ny(P_cc.ny)
    {
        rho.resize(nx, std::vector<double>(ny));
        rhovx.resize(nx, std::vector<double>(ny));
        rhovy.resize(nx, std::vector<double>(ny));
        rhovz.resize(nx, std::vector<double>(ny));
        Bx.resize(nx, std::vector<double>(ny));
        By.resize(nx, std::vector<double>(ny));
        Bz.resize(nx, std::vector<double>(ny));

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                rho[i][j] = P_cc.rho[i][j];
                rhovx[i][j] = P_cc.rho[i][j] * P_cc.vx[i][j];
                rhovy[i][j] = P_cc.rho[i][j] * P_cc.vy[i][j];
                rhovz[i][j] = P_cc.rho[i][j] * P_cc.vz[i][j];
                Bx[i][j] = P_cc.Bx[i][j];
                By[i][j] = P_cc.By[i][j];
                Bz[i][j] = P_cc.Bz[i][j];
            }
        }
        P = P_cc.P;
    }

    ~ConservativeVariablesCC() = default;

    void set(const ReconstructedValues& rv, int i, int j){
        rho[i][j] = rv.rho;
        rhovx[i][j] = rv.vx;
        rhovy[i][j] = rv.vy;
        rhovz[i][j] = rv.vz;
        Bx[i][j] = rv.Bx;
        By[i][j] = rv.By;
        Bz[i][j] = rv.Bz;
    }

    ReconstructedValues operator()(int i, int j) const {
        if (i < 0 || i >= nx || j < 0 || j >= ny)
            throw("Index out of range");

        return ReconstructedValues{rho[i][j], rhovx[i][j], rhovy[i][j], rhovz[i][j], Bx[i][j], By[i][j], Bz[i][j], P};
    }

    ConservativeVariablesCC operator*(double scalar) const {
        ConservativeVariablesCC result(this->nx, this->ny);

        for (int i = 0; i < this->nx; ++i) {
            for (int j = 0; j < this->ny; ++j) {
                result.rho[i][j] = this->rho[i][j] * scalar;
                result.rhovx[i][j] = this->rhovx[i][j] * scalar;
                result.rhovy[i][j] = this->rhovy[i][j] * scalar;
                result.rhovz[i][j] = this->rhovz[i][j] * scalar;
                result.Bx[i][j] = this->Bx[i][j] * scalar;
                result.By[i][j] = this->By[i][j] * scalar;
                result.Bz[i][j] = this->Bz[i][j] * scalar;
            }
        }
        return result;
    }

    ConservativeVariablesCC operator-(const ConservativeVariablesCC& other) const {
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

    ConservativeVariablesCC operator+(const ConservativeVariablesCC& other) const {
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

    ConservativeVariablesCC& operator=(const ConservativeVariablesCC& other) {
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
};

#endif // CONSERVATIVE_VARIABLES_CC_HPP_
