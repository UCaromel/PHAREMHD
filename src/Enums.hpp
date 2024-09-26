#ifndef ENUMS_HPP_
#define ENUMS_HPP_

enum BoundaryConditions { Periodic, ZeroGradient };

enum Reconstruction { Constant, Linear, WENO3, WENO5, WENO7 };
enum Slope { VanLeer, MinMod };
enum Riemann { Rusanov, HLL };
enum CTMethod { Average, Arithmetic, Contact, UCT_HLL, NoCT };
enum Integrator { EulerIntegrator, TVDRK2Integrator, TVDRK3Integrator };

enum OptionalPhysics {
  Hall,
  Res,
  Hyper,
  HallRes,
  HallHyper,
  ResHyper,
  HallResHyper,
  Off
};

enum dumpVariables { Primitive, Conservative, Both };

#endif
