import pyMHD as p
import initial_conditions as i

import os

result_dir = 'results/'
os.makedirs(result_dir, exist_ok=True)

P0cc = p.PrimitiveVariablesCC(i.nx, i.ny)
P0cc.init(i.rho, i.vx, i.vy, i.vz,i.Bx, i.By, i.Bz, i.P)

p.PhareMHD(P0cc, result_dir, i.order, i.nghost, 
           i.reconstruction, i.riemannsolver, i.constainedtransport, i.timeintegrator,
           i.Dx, i.Dy, i.FinalTime)