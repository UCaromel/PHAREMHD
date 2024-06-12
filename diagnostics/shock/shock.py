import sys
sys.path.append('/home/caromel/Documents/MHD_PHARE/pyMHD')

import numpy as np
import pyMHD as p
import os
import shutil

#############################################################################################################################################################################

nx = 800
ny = 1
Dx = 1
Dy = 1
Dt = 0.2
FinalTime = 80
order = 1
nghost = 2

boundaryconditions = p.BoundaryConditions.ZeroGradient

reconstruction = p.Reconstruction.Constant
slopelimiter = p.Slope.MinMod
riemannsolver = p.RiemannSolver.HLL
constainedtransport = p.CTMethod.Average
timeintegrator = p.Integrator.TVDRK2Integrator

dumpfrequency = 10
dumpvariables = p.dumpVariables.Primitive

##############################################################################################################################################################################

def rho_(x, y):
    return np.where(x<(nx*Dx/2), 1, 0.125)

def vx_(x, y):
    return 0.0

def vy_(x, y):
    return 0.0

def vz_(x, y):
    return 0.0

def Bx_(x, y):
    return 0.75

def By_(x, y):
    return np.where(x<(nx*Dx/2), 1, -1)

def Bz_(x, y):
    return 0.0

def P_(x, y):
    return np.where(x<(nx*Dx/2), 1, 0.1)

x = np.arange(nx) * Dx + 0.5 * Dx
y = np.arange(ny) * Dy + 0.5 * Dy

xx, yy = np.meshgrid(x, y, indexing = 'ij')

rho = np.full((nx, ny), rho_(xx, yy)).T
vx = np.full((nx, ny), vx_(xx, yy)).T
vy = np.full((nx, ny), vy_(xx, yy)).T
vz = np.full((nx, ny), vz_(xx, yy)).T
Bx = np.full((nx, ny), Bx_(xx, yy)).T
By = np.full((nx, ny), By_(xx, yy)).T
Bz = np.full((nx, ny), Bz_(xx, yy)).T
P = np.full((nx, ny), P_(xx, yy)).T

#############################################################################################################################################################################

result_dir = 'shockres/'
if os.path.exists(result_dir):
    shutil.rmtree(result_dir)

os.makedirs(result_dir, exist_ok=True)

P0cc = p.PrimitiveVariablesCC(nx, ny)
P0cc.init(rho, vx, vy, vz, Bx, By, Bz, P)

p.PhareMHD(P0cc, result_dir, order, nghost, 
           boundaryconditions, reconstruction, slopelimiter, riemannsolver, constainedtransport, timeintegrator,
           Dx, Dy, FinalTime, Dt, dumpvariables = dumpvariables)