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
        Dx = 0.01;
        Dy = 0.01;
        Dt = 0.00;
        FinalTime = 1;
        order = 1;
        nghost = 1;


        //double kx = (2*M_PI/10) * Dx + Dx/2.0;
        //double ky = (2*M_PI/10) * Dy + Dy/2.0;
        double k = (2.0*M_PI);
        double alpha = M_PI/4.0;
        double sinalpha = std::sin(alpha) * (2.0/std::sqrt(2));
        double cosalpha = std::cos(alpha) * (2.0/std::sqrt(2));

        rho.resize(ny, std::vector<double>(nx, 1.0));
        vx.resize(ny, std::vector<double>(nx, 0.0));
        vy.resize(ny, std::vector<double>(nx, 0.0));
        vz.resize(ny, std::vector<double>(nx, 0.0));
        Bx.resize(ny, std::vector<double>(nx, 1.0/(2.0*sinalpha)));
        By.resize(ny, std::vector<double>(nx, 1.0/(2.0*cosalpha)));
        Bz.resize(ny, std::vector<double>(nx, 0.0));
        P.resize(ny, std::vector<double>(nx, 0.1));

        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                vx[j][i] = UserFunction1(0.1, k*(0.5 + i)*cosalpha*Dx + k*(0.5 + j)*sinalpha*Dy) / 2.0 * sinalpha;
                vy[j][i] = UserFunction1(-0.1, k*(0.5 + i)*cosalpha*Dx + k*(0.5 + j)*sinalpha*Dy) / 2.0 * cosalpha;
            }
        }

        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                Bx[j][i] -= UserFunction1(0.1, k*(0.5 + i)*cosalpha*Dx + k*(0.5 + j)*sinalpha*Dy) / 2.0 * sinalpha;
                By[j][i] += UserFunction1(0.1, k*(0.5 + i)*cosalpha*Dx + k*(0.5 + j)*sinalpha*Dy) / 2.0 * cosalpha;
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