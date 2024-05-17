import pyMHD as p
import initial_conditions as i
import writting_utils as w

import numpy as np
import os

result_dir = 'results/'
os.makedirs(result_dir, exist_ok=True)

P0cc = p.PrimitiveVariablesCC(i.nx, i.ny)
P0cc.init(i.rho, i.vx, i.vy, i.vz,i.Bx, i.By, i.Bz, i.P)

P0 = p.InitialiseGhostCells(P0cc, i.nghost)

w.savePrimitiveVariables(P0, result_dir + 'Pcc_0.txt', i.nghost)

U0 = p.ConservativeVariablesCC(P0)
w.saveConcervativeVariables(U0, result_dir + "URK2_0.txt", i.nghost)
p.UpdateGhostCells(U0, i.nghost)

Un1.assign(U0)