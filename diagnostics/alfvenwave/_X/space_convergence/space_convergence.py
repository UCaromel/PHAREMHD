import sys
sys.path.append('/home/caromel/Documents/MHD_PHARE/pyMHD')

import numpy as np
import pyMHD as p
import os
import shutil

# Initial parameters
initial_nx = 50
initial_Dx = 0.02
final_time = 1
Dt = 5e-4
ny = 3
Dy = 0.01
order = 1
nghost = 2
reconstruction = p.Reconstruction.Linear
slopelimiter = p.Slope.VanLeer
riemannsolver = p.RiemannSolver.Rusanov
constainedtransport = p.CTMethod.Average
timeintegrator = p.Integrator.TVDRK2Integrator
dump_frequency = 10

###########################################################################################################################################################################

kx = 2 * np.pi

def rho_(x, y):
    return 1.0

def vx_(x, y):
    return 0.0

def vy_(x, y):
    return 1e-6 * np.cos(kx * x)

def vz_(x, y):
    return 0.0

def Bx_(x, y):
    return 1.0

def By_(x, y):
    return -1e-6 * np.cos(kx * x)

def Bz_(x, y):
    return 0.0

def P_(x, y):
    return 0.1


def initialize_variables(nx, ny, Dx, Dy):
    x = np.arange(nx) * Dx + 0.5 * Dx
    y = np.arange(ny) * Dy + 0.5 * Dy
    xx, yy = np.meshgrid(x, y, indexing='ij')
    rho = np.full((nx, ny), rho_(xx, yy)).T
    vx = np.full((nx, ny), vx_(xx, yy)).T
    vy = np.full((nx, ny), vy_(xx, yy)).T
    vz = np.full((nx, ny), vz_(xx, yy)).T
    Bx = np.full((nx, ny), Bx_(xx, yy)).T
    By = np.full((nx, ny), By_(xx, yy)).T
    Bz = np.full((nx, ny), Bz_(xx, yy)).T
    P = np.full((nx, ny), P_(xx, yy)).T
    return rho, vx, vy, vz, Bx, By, Bz, P

############################################################################################################################################################################
results_dir = 'space_results/'

if os.path.exists(results_dir):
    shutil.rmtree(results_dir)
os.makedirs(results_dir, exist_ok=True)

Dx = initial_Dx
nx = initial_nx
stepDx = 0

while Dx > initial_Dx / 32.0 and nx < 1600:
    rho, vx, vy, vz, Bx, By, Bz, P = initialize_variables(nx, ny, Dx, Dy)
    P0cc = p.PrimitiveVariablesCC(nx, ny)
    P0cc.init(rho, vx, vy, vz, Bx, By, Bz, P)

    current_results_dir = os.path.join(results_dir, f'{stepDx}_')
    
    p.PhareMHD(P0cc, current_results_dir, order, nghost, 
               reconstruction, slopelimiter, riemannsolver, constainedtransport, 
               timeintegrator, Dx, Dy, final_time, Dt, dump_frequency)
    
    stepDx += 1
    Dx /= 2.0
    nx *= 2
