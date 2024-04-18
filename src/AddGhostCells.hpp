#ifndef ADD_GHOST_CELLS_HPP_
#define ADD_GHOST_CELLS_HPP_

#include "ReconstructedValues.hpp"
#include "PrimitiveVariablesCC.hpp"
#include "ConservativeVariablesCC.hpp"


template<typename Variables>
Variables AddGhostCells(const Variables& V_cc, int nghost) {
    Variables withGhost(V_cc.nx + 2 * nghost, V_cc.ny + 2 * nghost);

    // Copy interior values
    for (int i = 0; i < V_cc.nx; ++i) {
        for (int j = 0; j < V_cc.ny; ++j) {
            withGhost.set(V_cc(i, j), i + nghost, j + nghost);
        }
    }

    if(nghost>V_cc.nx || nghost>V_cc.ny){
        throw("nghost too big for the domain");
    }

    if (nghost == 0) {
        return withGhost;
    } else {

        // For periodic boundary conditions

        // Add ghost cells for lines (excluding corners)
        for (int i = 0; i < V_cc.nx; i++) {
            for (int k = 1; k <= nghost; k++) {
                withGhost.set(V_cc(i, V_cc.ny - k), i + nghost, nghost - k); // Bottom side
                withGhost.set(V_cc(i, k - 1), i + nghost, withGhost.ny - 1 + k - nghost); // Top side
            }
        }

        // Add ghost cells for columns (excluding corners)
        for (int j = 0; j < V_cc.ny; j++) {
            for (int k = 1; k <= nghost; k++) {
                withGhost.set(V_cc(V_cc.nx - k, j), nghost - k, j + nghost); // Left side
                withGhost.set(V_cc(0, j), withGhost.nx - 1 + k - nghost, j + nghost); // Right side
            }
        }

        // Corners
        for (int k1 = 0; k1 < nghost; k1++) {
            for (int k2 = 0; k2 < nghost; k2++) {
                withGhost.set(0.5*(withGhost(k1 + V_cc.nx, k2) + withGhost(k1, k2 + V_cc.ny)), k1, k2); // Bottom-left
                withGhost.set(0.5*(withGhost(k1 + V_cc.nx, k2 + withGhost.ny - nghost) + withGhost(k1, k2 + nghost)), k1, k2 + withGhost.ny - nghost); // Top-left
                withGhost.set(0.5*(withGhost(k1 + nghost, k2) + withGhost(k1 + withGhost.nx - nghost, k2 + V_cc.ny)), k1 + withGhost.nx - nghost, k2); // Bottom-right
                withGhost.set(0.5*(withGhost(k1 + nghost, k2 + withGhost.ny - nghost) + withGhost(k1 + withGhost.nx - nghost, k2 + nghost)), k1 + withGhost.nx - nghost, k2 + withGhost.ny - nghost); // Top-right
            }
        }
        return withGhost;
    }
}



#endif //ADD_GHOST_CELLS_HPP_