#ifndef INTERFACE_HPP_
#define INTERFACE_HPP_

#include <cmath>

#include "ReconstructedValues.hpp"
#include "PrimitiveVariablesCC.hpp"
#include "EquationOfState.hpp"

enum struct Dir {
    X,
    Y
};

inline ReconstructedValues ComputeFluxVector(ReconstructedValues u, Dir dir) {
    ReconstructedValues flux;
    double GeneralisedPressure = u.P + 0.5*(u.Bx * u.Bx + u.By * u.By + u.Bz * u.Bz);

    if (dir == Dir::X) {
        flux.rho = u.rho * u.vx;
        flux.vx = u.rho * u.vx * u.vx + GeneralisedPressure - u.Bx * u.Bx;
        flux.vy = u.rho * u.vx * u.vy - u.Bx * u.By;
        flux.vz = u.rho * u.vx * u.vz - u.Bx * u.Bz;
        flux.Bx = 0.0;
        flux.By = u.By * u.vx - u.vy * u.Bx;
        flux.Bz = u.Bz * u.vx - u.vz * u.Bx;
        flux.P = (EosEtot(u) + GeneralisedPressure)*u.vx - u.Bx*(u.vx*u.Bx + u.vy*u.By + u.vz*u.Bz);
    } else if (dir == Dir::Y) {
        flux.rho = u.rho * u.vy;
        flux.vx = u.rho * u.vy * u.vx - u.By * u.Bx;
        flux.vy = u.rho * u.vy * u.vy + GeneralisedPressure - u.By * u.By;
        flux.vz = u.rho * u.vy * u.vz - u.By * u.Bz;
        flux.Bx = u.vx * u.By - u.vy * u.Bx;
        flux.By = 0.0;
        flux.Bz = u.Bz * u.vy - u.vz * u.By;
        flux.P = (EosEtot(u) + GeneralisedPressure)*u.vy - u.By*(u.vx*u.Bx + u.vy*u.By + u.vz*u.Bz);
    }
    return flux;
}

class Interface {
public:
    ReconstructedValues uL, uR, fL, fR;
    double SL, SR, Splus;

    // (i,j) interface index.
    Interface(const PrimitiveVariablesCC& P_cc /* Assuming ghost cells are added */, int i /* (0 to nx) + nghost */, int j /* (0 to ny) + nghost */, int order, int nghost, Dir dir) {
        if (order == 1) {
            if (dir == Dir::X) {
                uL = P_cc(i-1,j);
                uR = P_cc(i,j);
            } else if (dir == Dir::Y) {
                uL = P_cc(i,j-1);
                uR = P_cc(i,j);
            }
        }
        fL = ComputeFluxVector(uL, dir);
        fR = ComputeFluxVector(uR, dir);

        double gam = 5.0/3.0;

        double c0L = std::sqrt((gam*uL.P)/uL.rho); // Sound speeds
        double c0R = std::sqrt((gam*uR.P)/uR.rho); 

        double caxL = std::sqrt((uL.Bx*uL.Bx)/(uL.rho)); // Alfven speeds in x
        double caxR = std::sqrt((uR.Bx*uR.Bx)/(uR.rho));
        double cayL = std::sqrt((uL.By*uL.By)/(uL.rho)); // Alfven speeds in y
        double cayR = std::sqrt((uR.By*uR.By)/(uR.rho));

        double caL = std::sqrt((uL.Bx*uL.Bx + uL.By*uL.By + uL.Bz*uL.Bz)/(uL.rho)); // Alfven speeds
        double caR = std::sqrt((uR.Bx*uR.Bx + uR.By*uR.By + uR.Bz*uR.Bz)/(uR.rho));

        double cfastxL = std::sqrt((c0L*c0L + caL*caL)*0.5 + (std::sqrt((c0L*c0L + caL*caL)*(c0L*c0L + caL*caL) - 4*c0L*c0L*caxL*caxL))*0.5); // Fast magnetosonic speeds in x
        double cfastxR = std::sqrt((c0R*c0R + caR*caR)*0.5 + (std::sqrt((c0R*c0R + caR*caR)*(c0R*c0R + caR*caR) - 4*c0R*c0R*caxR*caxR))*0.5);
        double cfastyL = std::sqrt((c0L*c0L + caL*caL)*0.5 + (std::sqrt((c0L*c0L + caL*caL)*(c0L*c0L + caL*caL) - 4*c0L*c0L*cayL*cayL))*0.5); // Fast magnetosonic speeds in y
        double cfastyR = std::sqrt((c0R*c0R + caR*caR)*0.5 + (std::sqrt((c0R*c0R + caR*caR)*(c0R*c0R + caR*caR) - 4*c0R*c0R*cayR*cayR))*0.5);

        // Wave speeds
        if(dir == Dir::X){
            SL = std::min(uL.vx - cfastxL, uR.vx - cfastxR);
            SR = std::max(uL.vx + cfastxL, uR.vx + cfastxR);
            // For rusanov :
            Splus = std::max(std::abs(uL.vx) + cfastxL, std::abs(uR.vx) + cfastxR);
        }else if(dir == Dir::Y){
            SL = std::min(uL.vy - cfastyL, uR.vy - cfastyR);
            SR = std::max(uL.vy + cfastyL, uR.vy + cfastyR);
            //  For rusanov :
            Splus = std::max(std::abs(uL.vy) + cfastyL, std::abs(uR.vy) + cfastyR);
        }

        // Pass uL and uR in conservative form
        uL.vx = uL.rho * uL.vx;
        uL.vy = uL.rho * uL.vy;
        uL.vz = uL.rho * uL.vz;
        uL.P = EosEtot(uL);
        uR.vx = uR.rho * uR.vx;
        uR.vy = uR.rho * uR.vy;
        uR.vz = uR.rho * uR.vz;
        uR.P = EosEtot(uR);
    }
    ~Interface() = default;
};

#endif // INTERFACE_HPP_