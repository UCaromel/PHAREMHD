#ifndef INITIALISATION_HPP_
#define INITIALISATION_HPP_

#include <vector>
#include <cmath>

class Initialisation{
public:
    int nx;
    int ny;
    double Dx;
    double Dy;
    double Dt;
    double FinalTime;
    int order;
    int nghost;

    std::vector<std::vector<double>> rho;
    std::vector<std::vector<double>> vx;
    std::vector<std::vector<double>> vy;
    std::vector<std::vector<double>> vz;
    std::vector<std::vector<double>> Bx;
    std::vector<std::vector<double>> By;
    std::vector<std::vector<double>> Bz;
    double P;


    Initialisation() {
        nx = 100;
        ny = 100;
        Dx = 0.1;
        Dy = 0.1;
        Dt = 0.01;
        FinalTime = 1;
        order = 1;
        nghost = 1;
        rho.resize(nx, std::vector<double>(ny, (5.0/3.0)*(5.0/3.0)/4*M_PI));
        vx.resize(nx, std::vector<double>(ny, 1.0));
        vy.resize(nx, std::vector<double>(ny, 1.0));  // Corrected resize
        vz.resize(nx, std::vector<double>(ny, 1.0));
        Bx.resize(nx, std::vector<double>(ny, 1.0));  // Corrected resize
        By.resize(nx, std::vector<double>(ny, 1.0));
        Bz.resize(nx, std::vector<double>(ny, 1.0));  // Corrected resize
        //P = (5.0/3.0)/4*M_PI;
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                By[j][i] = UserFunction1(0.01, 2*M_PI, Dx*i);
            }
        }
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                vy[j][i] = UserFunction1(0.004, 2*M_PI, Dx*i);
            }
        }
    }
    ~Initialisation() = default;
private:
    double UserFunction1(double ampl, double k, double x){
        return ampl*cos(k*x);
    }
};


#endif //INITIALISATION_HPP_