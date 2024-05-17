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
    std::vector<std::vector<double>> P;


    Initialisation() {
        nx = 100;
        ny = 100;
        Dx = 0.1;
        Dy = 0.1;
        Dt = 0.00;
        FinalTime = 10;
        order = 1;
        nghost = 1;
        rho.resize(ny, std::vector<double>(nx, 1.0));
        vx.resize(ny, std::vector<double>(nx, 0.0));
        vy.resize(ny, std::vector<double>(nx, 0.0));
        vz.resize(ny, std::vector<double>(nx, 0.0));
        Bx.resize(ny, std::vector<double>(nx, 1.0));
        By.resize(ny, std::vector<double>(nx, 0.0));
        Bz.resize(ny, std::vector<double>(nx, 0.0));
        P.resize(ny, std::vector<double>(nx, 0.1));


        //double kx = (2*M_PI/10) * Dx + Dx/2.0;
        //double ky = (2*M_PI/10) * Dy + Dy/2.0;
        double k = (2*M_PI/10);
        double alpha = 45 * M_PI/180;

        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                vx[j][i] = UserFunction1(-1e-6, k*(0.5 + i)*std::cos(alpha)*Dx + k*(0.5 + j)*std::sin(alpha)*Dy) * std::sin(alpha);
                vy[j][i] = UserFunction1(-1e-6, k*(0.5 + i)*std::cos(alpha)*Dx + k*(0.5 + j)*std::sin(alpha)*Dy) * std::cos(alpha);
            }
        }

        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                Bx[j][i] += UserFunction1(1e-6, k*(0.5 + i)*std::cos(alpha)*Dx + k*(0.5 + j)*std::sin(alpha)*Dy) * std::sin(alpha);
                By[j][i] += UserFunction1(1e-6, k*(0.5 + i)*std::cos(alpha)*Dx + k*(0.5 + j)*std::sin(alpha)*Dy) * std::cos(alpha);
            }
        }
    }
    ~Initialisation() = default;
private:
    double UserFunction1(double ampl, double kx){
        return ampl*sin(kx);
    }
};


#endif //INITIALISATION_HPP_