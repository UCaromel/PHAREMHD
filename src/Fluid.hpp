#ifndef FLUID_HPP_
#define FLUID_HPP_

#include <array>

class Fluid
{
private:
    std::array<double,4> Pc; //primitive cell-centered
    std::array<double,4> Pf; //face-centered
    std::array<double,4> Pe; //edge-centered
    std::array<double,4> Cc; //conservative cell-centered
public:
    Fluid();
    ~Fluid();
};

#endif