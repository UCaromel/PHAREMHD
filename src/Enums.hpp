#ifndef ENUMS_HPP_
#define ENUMS_HPP_

enum Reconstruction { Constant, Linear };
enum Slope { VanLeer };
enum Riemann { Rusanov, HLL };
enum CTMethod { Average, UCT_HLL };
enum Integrator { EulerIntegrator, TVDRK2Integrator, TVDRK3Integrator };


#endif