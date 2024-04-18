#ifndef INTERFACE_HPP_
#define INTERFACE_HPP_

#include <cmath>

#include "ReconstructedValues.hpp"
#include "PrimitiveVariablesCC.hpp"
#include "AddGhostCells.hpp"

//might need a conversion to conserved variables here.

enum struct Dir {
    X,
    Y
};

ReconstructedValues ComputeFluxVector(ReconstructedValues u, Dir dir) {
    ReconstructedValues flux;
    if (dir == Dir::X) {
        flux.rho = u.rho * u.vx;
        flux.vx = u.rho * u.vx * u.vx + u.P - u.Bx * u.Bx;
        flux.vy = u.rho * u.vx * u.vy - u.Bx * u.By;
        flux.vz = u.rho * u.vx * u.vz - u.Bx * u.Bz;
        flux.Bx = 0.0;
        flux.By = u.By * u.vx - u.vy * u.Bx;
        flux.Bz = u.Bz * u.vx - u.vz * u.Bx;
        //flux.P = u.P;
    } else if (dir == Dir::Y) {
        flux.rho = u.rho * u.vy;
        flux.vx = u.rho * u.vy * u.vx - u.By * u.Bx;
        flux.vy = u.rho * u.vy * u.vy + u.P - u.By * u.By;
        flux.vz = u.rho * u.vy * u.vz - u.By * u.Bz;
        flux.Bx = u.vx * u.By - u.vy * u.Bx;
        flux.By = 0.0;
        flux.Bz = u.Bz * u.vy - u.vz * u.By;
        //flux.P = u.P;
    }
    return flux;
}

class Interface {
public:
    ReconstructedValues uL, uR, fL, fR;
    double SL, SR, Splus;

    // (i,j) interface index.
    Interface(PrimitiveVariablesCC& P_cc /* Assuming ghost cells are added */, int i /* 0 to nx */, int j /* 0 to ny */, int order, int nghost, Dir dir) {
        if (order == 1) {
            if (dir == Dir::X) {
                i += nghost;
                j += nghost;
                uL = P_cc(i-1,j);
                uR = P_cc(i,j);
            } else if (dir == Dir::Y) {
                i += nghost;
                j += nghost;
                uL = P_cc(i,j-1);
                uR = P_cc(i,j);
            }
        }
        fL = ComputeFluxVector(uL, dir);
        fR = ComputeFluxVector(uR, dir);

        double gamma = 5/3;
        double c0L = std::sqrt((gamma*uL.P)/uL.rho); // Sound speeds
        double c0R = std::sqrt((gamma*uR.P)/uR.rho); 
        double caxL = std::sqrt((uL.Bx*uL.Bx)/(4*M_PI*uL.rho)); // Alfven speeds in x
        double caxR = std::sqrt((uR.Bx*uR.Bx)/(4*M_PI*uR.rho));
        double caL = std::sqrt((uL.Bx*uL.Bx + uL.By*uL.By + uL.Bz*uL.Bz)/(4*M_PI*uL.rho)); // Alfven speeds
        double caR = std::sqrt((uR.Bx*uR.Bx + uR.By*uR.By + uR.Bz*uR.Bz)/(4*M_PI*uR.rho));
        double cfastL = std::sqrt((c0L*c0L + caL*caL)/2 + (std::sqrt((c0L*c0L + caL*caL)*(c0L*c0L + caL*caL) - 4*c0L*c0L*caxL*caxL))); // Fast magnetosonic speeds
        double cfastR = std::sqrt((c0R*c0R + caR*caR)/2 + (std::sqrt((c0R*c0R + caR*caR)*(c0R*c0R + caR*caR) - 4*c0R*c0R*caxR*caxR)));

        // Wave speeds
        if(dir == Dir::X){
            SL = std::min(uL.vx - cfastL, uR.vx - cfastR);
            SR = std::max(uL.vx + cfastL, uR.vx + cfastR);
            // For rusanov :
            Splus = std::max(std::abs(uL.vx) + cfastL, std::abs(uR.vx) + cfastR);
        }else if(dir == Dir::Y){
            SL = std::min(uL.vy - cfastL, uR.vy - cfastR);
            SR = std::max(uL.vy + cfastL, uR.vy + cfastR);
            //  For rusanov :
            Splus = std::max(std::abs(uL.vy) + cfastL, std::abs(uR.vy) + cfastR);
        }

        // Pass uL and uR in conservative form
        uL.vx = uL.rho * uL.vx;
        uL.vy = uL.rho * uL.vy;
        uL.vz = uL.rho * uL.vz;
        uR.vx = uR.rho * uR.vx;
        uR.vy = uR.rho * uR.vy;
        uR.vz = uR.rho * uR.vz;
    }
    ~Interface() = default;
};

#endif // INTERFACE_HPP_