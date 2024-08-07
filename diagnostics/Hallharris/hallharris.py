import sys
sys.path.append('/home/caromel/Documents/PHAREMHD/pyMHD')

import numpy as np
import pyMHD as p
import os
import shutil

#############################################################################################################################################################################

nx = 250
ny = 250
Dx = 0.1
Dy = 0.1
Dt = 0.0
FinalTime = 0.5
nghost = 2

boundaryconditions = p.BoundaryConditions.Periodic

reconstruction = p.Reconstruction.Linear
slopelimiter = p.Slope.MinMod
riemannsolver = p.RiemannSolver.Rusanov
constainedtransport = p.CTMethod.Average
timeintegrator = p.Integrator.TVDRK2Integrator

OptionalPhysics = p.OptionalPhysics.HallResHyper

dumpvariables = p.dumpVariables.Primitive
dumpfrequency = 1000

##############################################################################################################################################################################
Lx = nx*Dx
Ly = ny*Dy

def S(y, y0, l):
    return 0.5 * (1.0 + np.tanh((y - y0) / l))

def rho_(x, y):
    return 0.2 + 1.0 / np.cosh((y - Ly * 0.3) / 0.5) ** 2 + 1.0 / np.cosh((y - Ly * 0.7) / 0.5) ** 2

def vx_(x, y):
    return 0.0

def vy_(x, y):
    return 0.0

def vz_(x, y):
    return 0.0

def Bx_(x, y):
    w1 = 0.2
    w2 = 1.0
    x0 = x - 0.5 * Lx
    y1 = y - 0.3 * Ly
    y2 = y - 0.7 * Ly
    w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
    w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
    w5 = 2.0 * w1 / w2
    v1 = -1
    v2 = 1.0
    return (
        v1
        + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))
        + (-w5 * y1 * w3)
        + (+w5 * y2 * w4)
    )

def By_(x, y):
    w1 = 0.2
    w2 = 1.0
    x0 = x - 0.5 * Lx
    y1 = y - 0.3 * Ly
    y2 = y - 0.7 * Ly
    w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
    w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
    w5 = 2.0 * w1 / w2
    return (w5 * x0 * w3) + (-w5 * x0 * w4)

def Bz_(x, y):
    return 0.0


def P_(x, y):
    return 1.0 - (Bx_(x, y) ** 2)/2.0

x = np.arange(nx) * Dx + 0.5 * Dx
y = np.arange(ny) * Dy + 0.5 * Dy
xf = np.arange(nx+1) * Dx
yf = np.arange(ny+1) * Dy

xx, yy = np.meshgrid(x, y, indexing = 'ij')
xfx, yfx = np.meshgrid(xf, y, indexing = 'ij')
xfy, yfy = np.meshgrid(x, yf, indexing = 'ij')

rho = np.full((nx, ny), rho_(xx, yy)).T
vx = np.full((nx, ny), vx_(xx, yy)).T
vy = np.full((nx, ny), vy_(xx, yy)).T
vz = np.full((nx, ny), vz_(xx, yy)).T
Bxf = np.full((nx + 1, ny), Bx_(xfx, yfx)).T
Byf = np.full((nx, ny + 1), By_(xfy, yfy)).T
Bz = np.full((nx, ny), Bz_(xx, yy)).T
P = np.full((nx, ny), P_(xx, yy)).T

#############################################################################################################################################################################

result_dir = 'hallharrisresCT2/'
if os.path.exists(result_dir):
    shutil.rmtree(result_dir)

os.makedirs(result_dir, exist_ok=True)

P0cc = p.PrimitiveVariables(nx, ny)
P0cc.init(rho, vx, vy, vz, Bxf, Byf, Bz, P)

p.PhareMHD(P0cc, result_dir, nghost, 
           boundaryconditions, reconstruction, slopelimiter, riemannsolver, constainedtransport, timeintegrator,
           Dx, Dy, FinalTime, Dt, dumpvariables = dumpvariables, dumpfrequency = dumpfrequency, OptionalPhysics = OptionalPhysics)