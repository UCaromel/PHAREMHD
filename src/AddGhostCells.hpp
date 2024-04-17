#ifndef ADD_GHOST_CELLS_HPP_
#define ADD_GHOST_CELLS_HPP_

#include "ReconstructedValues.hpp"
#include "PrimitiveVariablesCC.hpp"

PrimitiveVariablesCC AddGhostCells(const PrimitiveVariablesCC& P_cc, int nghost) {
    PrimitiveVariablesCC withGhost(P_cc.nx + 2 * nghost, P_cc.ny + 2 * nghost);

    // Copy interior values
    for (int i = 0; i < P_cc.nx; ++i) {
        for (int j = 0; j < P_cc.ny; ++j) {
            withGhost.set(P_cc(i, j), i + nghost, j + nghost);
        }
    }

    if(nghost>P_cc.nx || nghost>P_cc.ny){
        throw("nghost too big for the domain");
    }

    if (nghost == 0) {
        return withGhost;
    } else {

        // For periodic boundary conditions

        // Add ghost cells for lines (excluding corners)
        for (int i = 0; i < P_cc.nx; ++i) {
            for (int k = 1; k <= nghost; ++k) {
                withGhost.set(P_cc(i, P_cc.ny - k), i + nghost, nghost - k); // Bottom side
                withGhost.set(P_cc(i, k - 1), i + nghost, withGhost.ny - 1 + k - nghost); // Top side
            }
        }

        // Add ghost cells for columns (excluding corners)
        for (int j = 0; j < P_cc.ny; ++j) {
            for (int k = 1; k < nghost; ++k) {
                withGhost.set(P_cc(P_cc.nx - k, j), nghost - k, j + nghost); // Left side
                withGhost.set(P_cc(0, j), withGhost.nx - 1 + k - nghost, j + nghost); // Right side
            }
        }

        // Corners are not useful for first order reconstruction

        return withGhost;
    }
}



#endif //ADD_GHOST_CELLS_HPP_