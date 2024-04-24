#include "EquationOfState.hpp"

const double gam = 3.0/5.0;

double EosEtot(const ReconstructedValues& rv){
    double Etot = rv.P/(gam - 1) + 0.5*(rv.rho*(rv.vx*rv.vx+ rv.vy*rv.vy + rv.vz*rv.vz)) + 0.5*(rv.Bx*rv.Bx + rv.By*rv.By+ rv.Bz*rv.Bz);
    return Etot;
}

double EosP(const ReconstructedValues& rv){
    double P = (gam - 1)*(rv.P - 0.5*((rv.vx*rv.vx + rv.vy*rv.vy + rv.vz*rv.vz)/rv.rho) - 0.5*(rv.Bx*rv.Bx + rv.By*rv.By + rv.Bz*rv.Bz));
    return P;
}