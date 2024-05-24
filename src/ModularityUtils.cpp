#include "ModularityUtils.hpp"

RiemannSolverFunction getRiemannSolver(Riemann rs){
    switch (rs)
    {
        case Rusanov: return &RusanovRiemannSolver;
        case HLL: return &HLLRiemannSolver;
        default: throw std::invalid_argument("Unknown riemann solver");
    }
}

IntegratorFunction getIntegrator(Integrator intg) {
    switch (intg) {
        case EulerIntegrator: return &Euler;
        case TVDRK2Integrator: return &TVDRK2;
        case TVDRK3Integrator: return &TVDRK3;
        default: throw std::invalid_argument("Unknown integrator");
    }
}