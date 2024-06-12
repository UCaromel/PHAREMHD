#ifndef ADD_GHOST_CELLS_HPP_
#define ADD_GHOST_CELLS_HPP_

#include <stdexcept>

#include "Enums.hpp"
#include "ReconstructedValues.hpp"
#include "PrimitiveVariablesCC.hpp"
#include "ConservativeVariablesCC.hpp"

template<typename Variables>
void UpdateGhostCells(Variables& V_cc, int nghost, BoundaryConditions bc) {
    if(bc == BoundaryConditions::Periodic){
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
        for (int k1 = 1; k1 <= nghost; k1++) {
            for (int k2 = 1; k2 <= nghost; k2++) {
                V_cc.set(V_cc(V_cc.nx - k1 - nghost, V_cc.ny - k2 - nghost), nghost - k1, nghost - k2); // Bottom-left
                V_cc.set(V_cc(V_cc.nx - k1 - nghost,  k2 - 1 + nghost), nghost - k1, k2 - 1 + V_cc.ny - nghost); // Top-left
                V_cc.set(V_cc(k1 - 1 + nghost, V_cc.ny - k2 - nghost), k1 - 1 + V_cc.nx - nghost, nghost - k2); // Bottom-right
                V_cc.set(V_cc(k1 - 1 + nghost, k2 - 1 + nghost), k1 - 1 + V_cc.nx - nghost, k2 - 1 + V_cc.ny - nghost); // Top-right
            }
        }
    } else if(bc == BoundaryConditions::ZeroGradient){
        // Update ghost cells for lines (excluding corners)
        for (int i = nghost; i < V_cc.nx - nghost; i++) {
            for (int k = 1; k <= nghost; k++) {
                V_cc.set(V_cc(i, nghost - k + 1), i, nghost - k); // Bottom side
                V_cc.set(V_cc(i, V_cc.ny - 2 + k - nghost), i, V_cc.ny - 1 + k - nghost); // Top side
            }
        }

        // Update ghost cells for columns (excluding corners)
        for (int j = nghost; j < V_cc.ny - nghost; j++) {
            for (int k = 1; k <= nghost; k++) {
                V_cc.set(V_cc(nghost - k + 1, j), nghost - k, j); // Left side
                V_cc.set(V_cc(V_cc.nx - 2 + k - nghost, j), V_cc.nx - 1 + k - nghost, j); // Right side
            }
        }

        // Corners (some more work for nghost > 1)
        for (int k1 = 1; k1 <= nghost; k1++) {
            for (int k2 = 1; k2 <= nghost; k2++) {
                V_cc.set(V_cc(nghost - k1 + 1, nghost - k2 + 1), nghost - k1, nghost - k2); // Bottom-left
                V_cc.set(V_cc(nghost - k1 + 1, k2 - 2 + V_cc.ny - nghost), nghost - k1, k2 - 1 + V_cc.ny - nghost); // Top-left
                V_cc.set(V_cc(k1 - 2 + V_cc.nx - nghost, nghost - k2), k1 - 1 + V_cc.nx - nghost, nghost - k2 + 1); // Bottom-right
                V_cc.set(V_cc(k1 - 2 + V_cc.nx - nghost, k2 - 2 + V_cc.ny - nghost), k1 - 1 + V_cc.nx - nghost, k2 - 1 + V_cc.ny - nghost); // Top-right
            }
        }
    } else {
        throw std::invalid_argument("Undefined boundary conditions");
    }
}

template<typename Variables>
Variables InitialiseGhostCells(const Variables& P_cc, int nghost, BoundaryConditions bc){
    Variables withGhost(P_cc.nx + 2 * nghost, P_cc.ny + 2 * nghost);

    // Copy interior values
    for (int i = 0; i < P_cc.nx; ++i) {
        for (int j = 0; j < P_cc.ny; ++j) {
            withGhost.set(P_cc(i, j), i + nghost, j + nghost);
        }
    }

    if (nghost == 0) {
        return withGhost;
    } else {
        UpdateGhostCells(withGhost, nghost, bc);
        return withGhost;
    }
}

#endif //ADD_GHOST_CELLS_HPP_