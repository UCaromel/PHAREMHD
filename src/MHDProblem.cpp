#include "Initialisation.hpp"
#include "ConservativeVariablesCC.hpp"
#include "PrimitiveVariablesCC.hpp"
#include "TimeIntegrator.hpp"


int main(){
    Initialisation I;

    PrimitiveVariablesCC P_cc(I.nx, I.ny);
    P_cc.init(I.rho, I.vx, I.vy, I.vz, I.Bx, I.By, I.Bz, I.P);

    for(double t=0.0; t<I.FinalTime; t += I.Dt){
        ConservativeVariablesCC U0(P_cc);
        ConservativeVariablesCC Un1 = TVDRK2(U0, I.Dx, I.Dy, I.Dt, I.order, I.nghost);
        P_cc = ConservativeVariablesCC(Un1);
    }
    
}