#ifndef ENUMS_HPP_
#define ENUMS_HPP_

enum Reconstruction { Constant };
enum Riemann { Rusanov, HLL };
enum CTMethod { Average };
enum Integrator { EulerIntegrator, TVDRK2Integrator, TVDRK3Integrator };

#endif