#ifndef SLOPE_LIMITER_HPP_
#define SLOPE_LIMITER_HPP_

#include <cmath>

#include "ReconstructedValues.hpp"

inline static ReconstructedValues VanLeerPhi(const ReconstructedValues& r) {
    return (r + r.rabs()) / (1.0 + r.rabs());
}

inline ReconstructedValues VanLeerSlope(const ReconstructedValues& Du12, const ReconstructedValues& Du_12){
    ReconstructedValues r;
    ReconstructedValues invr;

    r.rho = (Du12.rho == 0) ? 0 : Du_12.rho / Du12.rho;
    r.vx = (Du12.vx == 0) ? 0 : Du_12.vx / Du12.vx;
    r.vy = (Du12.vy == 0) ? 0 : Du_12.vy / Du12.vy;
    r.vz = (Du12.vz == 0) ? 0 : Du_12.vz / Du12.vz;
    r.Bx = (Du12.Bx == 0) ? 0 : Du_12.Bx / Du12.Bx;
    r.By = (Du12.By == 0) ? 0 : Du_12.By / Du12.By;
    r.Bz = (Du12.Bz == 0) ? 0 : Du_12.Bz / Du12.Bz;
    r.P = (Du12.P == 0) ? 0 : Du_12.P / Du12.P;

    invr.rho = (Du_12.rho == 0) ? 0 : Du12.rho / Du_12.rho;
    invr.vx = (Du_12.vx == 0) ? 0 : Du12.vx / Du_12.vx;
    invr.vy = (Du_12.vy == 0) ? 0 : Du12.vy / Du_12.vy;
    invr.vz = (Du_12.vz == 0) ? 0 : Du12.vz / Du_12.vz;
    invr.Bx = (Du_12.Bx == 0) ? 0 : Du12.Bx / Du_12.Bx;
    invr.By = (Du_12.By == 0) ? 0 : Du12.By / Du_12.By;
    invr.Bz = (Du_12.Bz == 0) ? 0 : Du12.Bz / Du_12.Bz;
    invr.P = (Du_12.P == 0) ? 0 : Du12.P / Du_12.P;

    ReconstructedValues Di = 0.5 * (1.0 + VanLeerPhi(r))*Du_12 + 0.5*(1.0 - VanLeerPhi(invr))*Du12;
    return Di;
}

#endif //SLOPE_LIMITER_HPP_