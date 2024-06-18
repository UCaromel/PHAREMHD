import sys
sys.path.append('/home/caromel/Documents/MHD_PHARE/pyMHD')

import numpy as np
import pyMHD as p
import os
import shutil

#############################################################################################################################################################################

nx = 100
ny = 100
Dx = 1/100
Dy = 1/100
Dt = 0.0
FinalTime = 0.15
order = 1
nghost = 2

boundaryconditions = p.BoundaryConditions.Periodic

reconstruction = p.Reconstruction.Linear
slopelimiter = p.Slope.VanLeer
riemannsolver = p.RiemannSolver.Rusanov
constainedtransport = p.CTMethod.Average
timeintegrator = p.Integrator.TVDRK2Integrator

##############################################################################################################################################################################
#omega = 20.0
B0 = 2.5/(np.sqrt(4*np.pi))#5./(np.sqrt(4*np.pi))
v0 = 1

r0 = 0.1
r1 = 0.115

# Domain [-0.5, 0.5]*[-0.5, 0.5]
def r(x, y):
    return np.sqrt((x-0.5)**2 + (y-0.5)**2)

def f(r):
    return (r1 - r)/(r1 - r0)

def rho_(x, y):
    r_ = r(x, y)
    f_ = f(r_)
    
    rho_values = np.where(r_ <= r0, 10.0, np.where(r_ < r1, 1.0 + 9.0 * f_, 1.0))
    return rho_values

def vx_(x, y):
    r_ = r(x, y)
    f_ = f(r_)
    
    vx_values = np.where(r_ <= r0, -v0*(y - 0.5)/r0, np.where(r_ < r1, -f_ * v0 * (y - 0.5)/ r_, 0.0))#np.where(r_ <= r0, -omega * y, np.where(r_ < r1, -f_ * omega * y * r0 / r_, 0.0))
    return vx_values

def vy_(x, y):
    r_ = r(x, y)
    f_ = f(r_)
    
    vy_values = np.where(r_ <= r0, v0*(x - 0.5)/r0, np.where(r_ < r1, f_ * v0 * (x - 0.5)/ r_, 0.0))#np.where(r_ <= r0, omega * x, np.where(r_ < r1, f_ * omega * x * r0 / r_, 0.0))
    return vy_values

def vz_(x, y):
    return 0.0

def Bx_(x, y):
    return B0

def By_(x, y):
    return 0.0

def Bz_(x, y):
    return 0.0 

def P_(x, y):
    return 1.0


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

result_dir = 'MHDrotorres/'
if os.path.exists(result_dir):
    shutil.rmtree(result_dir)

os.makedirs(result_dir, exist_ok=True)

P0cc = p.PrimitiveVariablesCC(nx, ny)
P0cc.init(rho, vx, vy, vz, Bx, By, Bz, P)

p.PhareMHD(P0cc, result_dir, order, nghost, 
           boundaryconditions, reconstruction, slopelimiter, riemannsolver, constainedtransport, timeintegrator,
           Dx, Dy, FinalTime, Dt)