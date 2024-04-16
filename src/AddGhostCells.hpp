#ifndef ADD_GHOST_CELLS_HPP_
#define ADD_GHOST_CELLS_HPP_

#include "ReconstructedValues.hpp"
#include "PrimitiveVariablesCC.hpp"

PrimitiveVariablesCC AddGhostCells(const PrimitiveVariablesCC& P_cc, int nghost){
    PrimitiveVariablesCC withGhost(P_cc.nx+2*nghost, P_cc.ny+2*nghost);

    //interior cells (non ghost)
    for(int i=0; i<P_cc.ny; i++){
        for(int j=0; j<P_cc.nx; j++){
            withGhost.set(P_cc[j],j+withGhost.nx+1+2*nghost*i);
        }
    }

    return withGhost;
}


#endif //ADD_GHOST_CELLS_HPP_

/*  //the following only works for one ghost cell

    //bottom ghost cells 
    for(int i=1; i<=P_cc.nx; i++){
        withGhost.set(P_cc[i+P_cc.nx*(P_cc.ny-1)-1],i);
    }

    //top ghost cells
    for(int i=withGhost.nx*withGhost.ny-P_cc.nx-1; i<=withGhost.nx*withGhost.ny-1;i++){
        withGhost.set(P_cc[i-(withGhost.nx*withGhost.ny-P_cc.nx-1)],i);
    }
*/