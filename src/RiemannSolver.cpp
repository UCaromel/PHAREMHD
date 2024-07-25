#include "RiemannSolver.hpp"

ReconstructedValues RusanovRiemannSolver(const Interface& inter){
    ReconstructedValues Frus;
    Frus = 0.5 * (inter.fL + inter.fR) - 0.5 * inter.Splus * (inter.uR - inter.uL);

    if (inter.OP == OptionalPhysics::HallResHyper){
        Frus.Bx = 0.5 * (inter.fL.Bx + inter.fR.Bx) - 0.5 * inter.Splusb * (inter.uR.Bx - inter.uL.Bx);
        Frus.By = 0.5 * (inter.fL.By + inter.fR.By) - 0.5 * inter.Splusb * (inter.uR.By - inter.uL.By);
        Frus.Bz = 0.5 * (inter.fL.Bz + inter.fR.Bz) - 0.5 * inter.Splusb * (inter.uR.Bz - inter.uL.Bz);
        Frus.P = 0.5 * (inter.fL.P + inter.fR.P) - 0.5 * inter.Splusb * (inter.uR.P - inter.uL.P);
    }

    return Frus;
}

ReconstructedValues HLLRiemannSolver(const Interface& inter){
    ReconstructedValues fhll;
    if(inter.SL > 0){
        fhll = inter.fL;
    }
    else if (inter.SR < 0){
        fhll = inter.fR;
    }
    else{
        fhll = (inter.SR * inter.fL - inter.SL * inter.fR + inter.SL * inter.SR * (inter.uR - inter.uL)) / (inter.SR - inter.SL);
    }


    if (inter.OP == OptionalPhysics::HallResHyper){
        if(inter.SLb > 0){
            fhll.Bx = inter.fL.Bx;
            fhll.By = inter.fL.By;
            fhll.Bz = inter.fL.Bz;
            fhll.P = inter.fL.P;
        }
        else if (inter.SRb < 0) {
            fhll.Bx = inter.fR.Bx;
            fhll.By = inter.fR.By;
            fhll.Bz = inter.fR.Bz;
            fhll.P = inter.fR.P;
        }
        else {
            fhll.Bx = (inter.SRb * inter.fL.Bx - inter.SLb * inter.fR.Bx + inter.SLb * inter.SRb * (inter.uR.Bx - inter.uL.Bx)) / (inter.SRb - inter.SLb);
            fhll.By = (inter.SRb * inter.fL.By - inter.SLb * inter.fR.By + inter.SLb * inter.SRb * (inter.uR.By - inter.uL.By)) / (inter.SRb - inter.SLb);
            fhll.Bz = (inter.SRb * inter.fL.Bz - inter.SLb * inter.fR.Bz + inter.SLb * inter.SRb * (inter.uR.Bz - inter.uL.Bz)) / (inter.SRb - inter.SLb);
            fhll.P = (inter.SRb * inter.fL.P - inter.SLb * inter.fR.P + inter.SLb * inter.SRb * (inter.uR.P - inter.uL.P)) / (inter.SRb - inter.SLb);
        }
    }

    return fhll;
}