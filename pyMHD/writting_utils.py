import numpy as np

def savePrimitiveVariables(P_cc, filename, nghost):
    # Extract the region without ghost cells and flatten each component
    rho = P_cc.rho[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    vx = P_cc.vx[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    vy = P_cc.vy[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    vz = P_cc.vz[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    Bx = P_cc.Bx[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    By = P_cc.By[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    Bz = P_cc.Bz[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    P = P_cc.P[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    
    # Stack columns
    data = np.column_stack((rho, vx, vy, vz, Bx, By, Bz, P))
    
    # Save to file
    np.savetxt(filename, data, fmt='%15.10f')
    print("Wrote file:", filename)

def saveConcervativeVariables(P_cc, filename, nghost):
    # Extract the region without ghost cells and flatten each component
    rho = P_cc.rho[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    rhovx = P_cc.rhovx[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    rhovy = P_cc.rhovy[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    rhovz = P_cc.rhovz[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    Bx = P_cc.Bx[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    By = P_cc.By[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    Bz = P_cc.Bz[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    Etot = P_cc.Etot[nghost:(P_cc.ny - nghost), nghost:(P_cc.nx - nghost)].flatten()
    
    # Stack columns
    data = np.column_stack((rho, rhovx, rhovy, rhovz, Bx, By, Bz, Etot))
    
    # Save to file
    np.savetxt(filename, data, fmt='%15.10f')
    print("Wrote file:", filename)


def saveVectorToFile(values, filename):
    with open(filename, 'w') as outFile:
        if isinstance(values[0], list):
            data = np.array(values)
            np.savetxt(outFile, data, fmt='%15.10f')
        else:
            data = np.array(values)
            np.savetxt(outFile, data.reshape(-1, 1), fmt='%15.10f')
        print("Wrote file:", filename)

def formatTime(time):
    return f"{time:.10f}".replace('.', '_')
