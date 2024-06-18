#include "ModularityUtils.hpp"

SLFunction getSlopeLimiter(Slope sl){
    switch (sl)
    {
        case VanLeer: return &VanLeerSlope;
        case MinMod: return &MinModSlope;
        default: throw std::invalid_argument("Unknown slope limiter");
    }
}

RiemannSolverFunction getRiemannSolver(Riemann rs){
    switch (rs)
    {
        case Rusanov: return &RusanovRiemannSolver;
        case HLL: return &HLLRiemannSolver;
        default: throw std::invalid_argument("Unknown riemann solver");
    }
}

CTFunction getCT(CTMethod ct){
    switch (ct){
        case Average: return &ConstrainedTransportAverage;
        case Contact: return &ConstrainedTransportContact;
        case UCT_HLL: return &UCTHLL;
        default: throw std::invalid_argument("Unknown constrained transport");
    }
}

IntegratorFunction getIntegrator(Integrator intg){
    switch (intg) {
        case EulerIntegrator: return &Euler;
        case TVDRK2Integrator: return &TVDRK2;
        case TVDRK3Integrator: return &TVDRK3;
        default: throw std::invalid_argument("Unknown integrator");
    }
}