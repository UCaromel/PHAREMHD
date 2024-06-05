import sys
sys.path.append('/home/caromel/Documents/MHD_PHARE/pyMHD')

import numpy as np
import pyMHD as p
import os
import shutil

#############################################################################################################################################################################

nx = 200
ny = 100
Dx = 0.01
Dy = 0.01
Dt = 0.0
FinalTime = 1
order = 1
nghost = 2

reconstruction = p.Reconstruction.Linear
slopelimiter = p.Slope.VanLeer
riemannsolver = p.RiemannSolver.HLL
constainedtransport = p.CTMethod.UCT_HLL
timeintegrator = p.Integrator.TVDRK2Integrator

##############################################################################################################################################################################
lx=nx*Dx
ly=ny*Dy

alpha = np.arctan(0.5)
cosalpha = np.cos(alpha)
sinalpha = np.sin(alpha)
kx = 2*np.pi/(cosalpha*lx)
ky = 2*np.pi/(sinalpha*ly)

def rho_(x, y):
    return 1.0

def vx_(x, y):
    return (-1e-4 * np.sin(kx * x * cosalpha + ky * y * sinalpha)) / (2. * sinalpha)

def vy_(x, y):
    return (1e-4 * np.sin(kx * x * cosalpha + ky * y * sinalpha)) / (2. * cosalpha)

def vz_(x, y):
    return 0.0

def Bx_(x, y):
    return (1 - 1e-4 * np.sin(kx * x * cosalpha + ky * y * sinalpha)) / (2. * sinalpha)

def By_(x, y):
    return (1 + 1e-4 * np.sin(kx * x * cosalpha + ky * y * sinalpha)) / (2. * cosalpha)

def Bz_(x, y):
    return 0.0

def P_(x, y):
    return 0.1

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

result_dir = '2Dwave/'
if os.path.exists(result_dir):
    shutil.rmtree(result_dir)

os.makedirs(result_dir, exist_ok=True)

P0cc = p.PrimitiveVariablesCC(nx, ny)
P0cc.init(rho, vx, vy, vz, Bx, By, Bz, P)

p.PhareMHD(P0cc, result_dir, order, nghost, 
           reconstruction, slopelimiter, riemannsolver, constainedtransport, timeintegrator,
           Dx, Dy, FinalTime)