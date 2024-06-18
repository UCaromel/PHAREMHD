#ifndef ENUMS_HPP_
#define ENUMS_HPP_

enum BoundaryConditions { Periodic, ZeroGradient };

enum Reconstruction { Constant, Linear };
enum Slope { VanLeer, MinMod };
enum Riemann { Rusanov, HLL };
enum CTMethod { Average, Contact, UCT_HLL };
enum Integrator { EulerIntegrator, TVDRK2Integrator, TVDRK3Integrator };

enum dumpVariables { Primitive, Conservative, Both };

#endif