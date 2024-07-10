#ifndef ADD_GHOST_CELLS_HPP_
#define ADD_GHOST_CELLS_HPP_

#include <stdexcept>

#include "Enums.hpp"
#include "ReconstructedValues.hpp"
#include "PrimitiveVariables.hpp"
#include "ConservativeVariables.hpp"

template<typename Variables>
void UpdateGhostCells(Variables& Vn, int nghost, BoundaryConditions bc) {
    if(bc == BoundaryConditions::Periodic){
        // Update ghost cells for lines (excluding corners)
        for (int i = nghost; i < Vn.nx - nghost; i++) {
            for (int k = 1; k <= nghost; k++) {
                Vn.set(Vn(i, Vn.ny - k - nghost), i, nghost - k); // Bottom side
                Vn.set(Vn(i, k - 1 + nghost), i, Vn.ny - 1 + k - nghost); // Top side
            }
        }

        // Update ghost cells for columns (excluding corners)
        for (int j = nghost; j < Vn.ny - nghost; j++) {
            for (int k = 1; k <= nghost; k++) {
                Vn.set(Vn(Vn.nx - k - nghost, j), nghost - k, j); // Left side
                Vn.set(Vn(k - 1 + nghost, j), Vn.nx - 1 + k - nghost, j); // Right side
            }
        }

        // Corners
        for (int k1 = 1; k1 <= nghost; k1++) {
            for (int k2 = 1; k2 <= nghost; k2++) {
                Vn.set(Vn(Vn.nx - k1 - nghost, Vn.ny - k2 - nghost), nghost - k1, nghost - k2); // Bottom-left
                Vn.set(Vn(Vn.nx - k1 - nghost,  k2 - 1 + nghost), nghost - k1, k2 - 1 + Vn.ny - nghost); // Top-left
                Vn.set(Vn(k1 - 1 + nghost, Vn.ny - k2 - nghost), k1 - 1 + Vn.nx - nghost, nghost - k2); // Bottom-right
                Vn.set(Vn(k1 - 1 + nghost, k2 - 1 + nghost), k1 - 1 + Vn.nx - nghost, k2 - 1 + Vn.ny - nghost); // Top-right
            }
        }

        // Face centered (only one ghost cell)
        for (int i = nghost; i < Vn.nx + 1 - nghost; i++) {
            Vn.Bxf[nghost - 1][i] = Vn.Bxf[Vn.ny - nghost - 1][i]; // Bottom side
            Vn.Bxf[Vn.ny - nghost][i] = Vn.Bxf[nghost][i]; // Top side
        }

        for (int j = nghost; j < Vn.ny + 1 - nghost; j++) {
            Vn.Byf[j][nghost - 1] = Vn.Byf[j][Vn.nx - nghost - 1]; // Left side
            Vn.Byf[j][Vn.nx - nghost] = Vn.Byf[j][nghost]; // Right side
        }

    } else if(bc == BoundaryConditions::ZeroGradient){
        // Update ghost cells for lines (excluding corners)
        for (int i = nghost; i < Vn.nx - nghost; i++) {
            for (int k = 1; k <= nghost; k++) {
                Vn.set(Vn(i, nghost - k + 1), i, nghost - k); // Bottom side
                Vn.set(Vn(i, Vn.ny - 2 + k - nghost), i, Vn.ny - 1 + k - nghost); // Top side
            }
        }

        // Update ghost cells for columns (excluding corners)
        for (int j = nghost; j < Vn.ny - nghost; j++) {
            for (int k = 1; k <= nghost; k++) {
                Vn.set(Vn(nghost - k + 1, j), nghost - k, j); // Left side
                Vn.set(Vn(Vn.nx - 2 + k - nghost, j), Vn.nx - 1 + k - nghost, j); // Right side
            }
        }

        // Corners
        for (int k1 = 1; k1 <= nghost; k1++) {
            for (int k2 = 1; k2 <= nghost; k2++) {
                Vn.set(Vn(nghost - k1 + 1, nghost - k2 + 1), nghost - k1, nghost - k2); // Bottom-left
                Vn.set(Vn(nghost - k1 + 1, k2 - 2 + Vn.ny - nghost), nghost - k1, k2 - 1 + Vn.ny - nghost); // Top-left
                Vn.set(Vn(k1 - 2 + Vn.nx - nghost, nghost - k2 + 1), k1 - 1 + Vn.nx - nghost, nghost - k2); // Bottom-right
                Vn.set(Vn(k1 - 2 + Vn.nx - nghost, k2 - 2 + Vn.ny - nghost), k1 - 1 + Vn.nx - nghost, k2 - 1 + Vn.ny - nghost); // Top-right
            }
        }

        // Face centered (only one ghost cell)
        for (int i = nghost; i < Vn.nx + 1 - nghost; i++) {
            Vn.Bxf[nghost - 1][i] = Vn.Bxf[nghost][i]; // Bottom side
            Vn.Bxf[Vn.ny - nghost][i] = Vn.Bxf[Vn.ny - nghost - 1][i]; // Top side
        }

        for (int j = nghost; j < Vn.ny + 1 - nghost; j++) {
            Vn.Byf[j][nghost - 1] = Vn.Byf[j][nghost]; // Left side
            Vn.Byf[j][Vn.nx - nghost] = Vn.Byf[j][Vn.nx - nghost - 1]; // Right side
        }

    } else {
        throw std::invalid_argument("Undefined boundary conditions");
    }
}

template<typename Variables>
Variables InitialiseGhostCells(const Variables& Vn, int nghost, BoundaryConditions bc){
    Variables withGhost(Vn.nx + 2 * nghost, Vn.ny + 2 * nghost);

    // Copy interior values
    for (int i = 0; i < Vn.nx; ++i) {
        for (int j = 0; j < Vn.ny; ++j) {
            withGhost.set(Vn(i, j), i + nghost, j + nghost);
        }
    }

    if (nghost == 0) {
        return withGhost;
    }

    // Face centered (only one ghost cell)
    for (int j = 0; j < Vn.ny; ++j) {
        for (int i = 0; i <= Vn.nx; ++i) {
            withGhost.Bxf[j + nghost][i + nghost] = Vn.Bxf[j][i];
        }
    }
    for (int j = 0; j <= Vn.ny; ++j) {
        for (int i = 0; i < Vn.nx; ++i) {
            withGhost.Byf[j + nghost][i + nghost] = Vn.Byf[j][i];
        }
    }

    UpdateGhostCells(withGhost, nghost, bc);
    return withGhost;
}

#endif //ADD_GHOST_CELLS_HPP_