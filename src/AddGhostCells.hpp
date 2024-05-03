#ifndef ADD_GHOST_CELLS_HPP_
#define ADD_GHOST_CELLS_HPP_

#include "ReconstructedValues.hpp"
#include "PrimitiveVariablesCC.hpp"
#include "ConservativeVariablesCC.hpp"

template<typename Variables>
void UpdateGhostCells(Variables& V_cc, int nghost) {
    // For periodic boundary conditions

    // Update ghost cells for lines (excluding corners)
    for (int i = nghost; i < V_cc.nx - nghost; i++) {
        for (int k = 1; k <= nghost; k++) {
            V_cc.set(V_cc(i, V_cc.ny - k - nghost), i, nghost - k); // Bottom side
            V_cc.set(V_cc(i, k - 1 + nghost), i, V_cc.ny - 1 + k - nghost); // Top side
        }
    }

    // Update ghost cells for columns (excluding corners)
    for (int j = nghost; j < V_cc.ny - nghost; j++) {
        for (int k = 1; k <= nghost; k++) {
            V_cc.set(V_cc(V_cc.nx - k - nghost, j), nghost - k, j); // Left side
            V_cc.set(V_cc(k - 1 + nghost, j), V_cc.nx - 1 + k - nghost, j); // Right side
        }
    }

    // Corners (some more work for nghost > 1)
    for (int k1 = 0; k1 < nghost; k1++) {
        for (int k2 = 0; k2 < nghost; k2++) {
            V_cc.set(0.5*(V_cc(k1 + V_cc.nx - 2*nghost, k2) + V_cc(k1, k2 + V_cc.ny - 2*nghost)), k1, k2); // Bottom-left
            V_cc.set(0.5*(V_cc(k1 + V_cc.nx - 2*nghost, k2 + V_cc.ny - nghost) + V_cc(k1, k2 + nghost)), k1, k2 + V_cc.ny - nghost); // Top-left
            V_cc.set(0.5*(V_cc(k1 + nghost, k2) + V_cc(k1 + V_cc.nx - nghost, k2 + V_cc.ny - 2*nghost)), k1 + V_cc.nx - nghost, k2); // Bottom-right
            V_cc.set(0.5*(V_cc(k1 + nghost, k2 + V_cc.ny - nghost) + V_cc(k1 + V_cc.nx - nghost, k2 + nghost)), k1 + V_cc.nx - nghost, k2 + V_cc.ny - nghost); // Top-right
        }
    }
/*
    for(int i=1; i < V_cc.nx - 1; ++i){
        V_cc.set(V_cc(i, V_cc.ny - 2), i, 0);
        V_cc.set(V_cc(i, 1), i, V_cc.ny - 1);
    }
    for(int j=1; j < V_cc.ny - 1; ++j){
        V_cc.set(V_cc(V_cc.nx - 2, j), 0, j);
        V_cc.set(V_cc(1, j), V_cc.nx - 1, j);
    }
    V_cc.set(V_cc(V_cc.nx-2,V_cc.ny-2),0,0);
    V_cc.set(V_cc(V_cc.nx-2,1),0, V_cc.ny-1);
    V_cc.set(V_cc(1,V_cc.ny-2),V_cc.nx-1,0);
    V_cc.set(V_cc(1,1),V_cc.nx-1,V_cc.ny-1);
*/
}

template<typename Variables>
Variables InitialiseGhostCells(const Variables& P_cc, int nghost){
    Variables withGhost(P_cc.nx + 2 * nghost, P_cc.ny + 2 * nghost);

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
        UpdateGhostCells(withGhost, nghost);
        return withGhost;
    }
}

#endif //ADD_GHOST_CELLS_HPP_