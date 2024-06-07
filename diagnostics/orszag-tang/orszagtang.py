import sys
sys.path.append('/home/caromel/Documents/MHD_PHARE/pyMHD')

import numpy as np
import pyMHD as p
import os
import shutil

#############################################################################################################################################################################

nx = 128
ny = 128
Dx = 1/128
Dy = 1/128
Dt = 0.0
FinalTime = 0.5
order = 1
nghost = 2

reconstruction = p.Reconstruction.Linear
slopelimiter = p.Slope.VanLeer
riemannsolver = p.RiemannSolver.Rusanov
constainedtransport = p.CTMethod.Average
timeintegrator = p.Integrator.TVDRK2Integrator

##############################################################################################################################################################################
lx=nx*Dx
ly=ny*Dy
B0 = 1./(np.sqrt(4.*np.pi))

def rho_(x, y):
    return 25./(36.*np.pi)

def vx_(x, y):
    return -np.sin(2.*np.pi*y)

def vy_(x, y):
    return np.sin(2.*np.pi*x)

def vz_(x, y):
    return 0.0

def Bx_(x, y):
    return -B0*np.sin(2.*np.pi*y)#0.0

def By_(x, y):
    return B0*np.sin(4.*np.pi*x)#0.0

def Bz_(x, y):
    return 0.0 #-np.sin(x) + np.sin(y)

def P_(x, y):
    return 5./(12.*np.pi)

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

result_dir = 'orszagtangres/'
if os.path.exists(result_dir):
    shutil.rmtree(result_dir)

os.makedirs(result_dir, exist_ok=True)

P0cc = p.PrimitiveVariablesCC(nx, ny)
P0cc.init(rho, vx, vy, vz, Bx, By, Bz, P)

p.PhareMHD(P0cc, result_dir, order, nghost, 
           reconstruction, slopelimiter, riemannsolver, constainedtransport, timeintegrator,
           Dx, Dy, FinalTime)