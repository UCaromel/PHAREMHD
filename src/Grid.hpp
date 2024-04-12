#ifndef GRID_HPP_
#define GRID_HPP_

#include <array>
#include <string>

#include "MHDProblem.hpp"

class Grid
{
public:
    std::array<std::array<double,1>,2> x;  //cell centers
    std::array<std::array<double,1>,2> xr;  //right interface
    std::array<std::array<double,1>,2> xl;  //left interface
    std::array<std::array<double,1>,2> dx;  //cell widths

    std::array<double,2> xbeg;  //beginning of the grid
    std::array<double,2> xend;  //end of the grid

    std::array<int,2> np_tot;  //total number of cells
    std::array<int,2> np_int;  //total number of cells without ghost cells
    std::array<int,2> nghost;  //total number of ghost cells

    std::array<BoundaryType,2> rbound;  //right boundary condition
    std::array<BoundaryType,2> lbound;  //left boundaty condition

    Grid(int nx, int ny, std::array<BoundaryType,2> r, std::array<BoundaryType,2> l);
};


#endif