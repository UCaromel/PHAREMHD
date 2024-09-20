import sys
sys.path.append('/home/caromel/Documents/PHAREMHD/pyMHD')

import numpy as np
import pyMHD as p
import os
import shutil

#############################################################################################################################################################################

nx = 128
ny = 1
Dx = 0.05
Dy = 1
Dt = 0.0006

nghost = 2

boundaryconditions = p.BoundaryConditions.Periodic

reconstruction = p.Reconstruction.Linear

slopelimiter = p.Slope.VanLeer
riemannsolver = p.RiemannSolver.Rusanov
constainedtransport = p.CTMethod.NoCT
timeintegrator = p.Integrator.EulerIntegrator

consts = p.Consts(sigmaCFL = 0.8, gam = 5/3, eta = 0.0, nu = 0.000)
physics = p.OptionalPhysics.HallResHyper

dumpvariables = p.dumpVariables.Primitive
dumpfrequency = 1

##############################################################################################################################################################################
lx=nx*Dx
k=2*np.pi/lx

m = 4

kt=2*np.pi/lx * m
w = (kt**2 /2) *(np.sqrt(1+4/kt**2) + 1)
FinalTime = 2*np.pi / w

np.random.seed(0)

modes = [4]
phases = np.random.rand(len(modes))

def rho_(x, y):
    return 1.0

def vx_(x, y):
    return 0.0

def vy_(x, y):
    ret = np.zeros((x.shape[0], y.shape[1]))
    
    for m,phi in zip(modes, phases):
        ret[:,:] += -np.cos(k*x*m + phi)*1e-2*k
    return ret

def vz_(x, y):
    ret = np.zeros((x.shape[0], y.shape[1]))
    for m,phi in zip(modes, phases):
        ret[:,:] += np.sin(k*x*m + phi)*1e-2*k
    return ret
def Bx_(x, y):
    return 1.0

def By_(x, y):
    ret = np.zeros((x.shape[0], y.shape[1]))
    for m,phi in zip(modes, phases):
        ret[:,:] += np.cos(k*x*m + phi)*1e-2
    return ret

def Bz_(x, y):
    ret = np.zeros((x.shape[0], y.shape[1]))
    for m,phi in zip(modes, phases):
        ret[:,:] += -np.sin(k*x*m + phi)*1e-2
    return ret

def P_(x, y):
    return 1.0

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

result_dir = 'whislerwaveres_singlemode/'
if os.path.exists(result_dir):
    shutil.rmtree(result_dir)

os.makedirs(result_dir, exist_ok=True)

P0cc = p.PrimitiveVariables(nx, ny)
P0cc.init(rho, vx, vy, vz, Bxf, Byf, Bz, P)

p.PhareMHD(P0cc, result_dir, nghost, 
           boundaryconditions, reconstruction, slopelimiter, riemannsolver, constainedtransport, timeintegrator,
           Dx, Dy, FinalTime, Dt = Dt, dumpvariables = dumpvariables, Consts = consts, OptionalPhysics = physics, dumpfrequency = dumpfrequency)