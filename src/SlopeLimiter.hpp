#ifndef SLOPE_LIMITER_HPP_
#define SLOPE_LIMITER_HPP_

#include <cmath>

#include "ReconstructedValues.hpp"

inline ReconstructedValues VanLeerSlope(const ReconstructedValues& Du12, const ReconstructedValues& Du_12){
    ReconstructedValues Di;

    Di.rho = Du12.rho * Du_12.rho > 0.0 ? 2.0 * Du12.rho * Du_12.rho / (Du12.rho + Du_12.rho) : 0.0;
    Di.vx = Du12.vx * Du_12.vx > 0.0 ? 2.0 * Du12.vx * Du_12.vx / (Du12.vx + Du_12.vx) : 0.0;
    Di.vy = Du12.vy * Du_12.vy > 0.0 ? 2.0 * Du12.vy * Du_12.vy / (Du12.vy + Du_12.vy) : 0.0;
    Di.vz = Du12.vz * Du_12.vz > 0.0 ? 2.0 * Du12.vz * Du_12.vz / (Du12.vz + Du_12.vz) : 0.0;
    Di.Bx = Du12.Bx * Du_12.Bx > 0.0 ? 2.0 * Du12.Bx * Du_12.Bx / (Du12.Bx + Du_12.Bx) : 0.0;
    Di.By = Du12.By * Du_12.By > 0.0 ? 2.0 * Du12.By * Du_12.By / (Du12.By + Du_12.By) : 0.0;
    Di.Bz = Du12.Bz * Du_12.Bz > 0.0 ? 2.0 * Du12.Bz * Du_12.Bz / (Du12.Bz + Du_12.Bz) : 0.0;
    Di.P = Du12.P * Du_12.P > 0.0 ? 2.0 * Du12.P * Du_12.P / (Du12.P + Du_12.P) : 0.0;
    return Di;
}

#endif //SLOPE_LIMITER_HPP_