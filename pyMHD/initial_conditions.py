import numpy as np

nx = 1000
ny = 10
Dx = 0.01
Dy = 0.01
Dt = 0.008
FinalTime = 2
order = 1
nghost = 1

def rho_(x, y):
    return 1.0

def vx_(x, y):
    return 0.0

def vy_(x, y):
    return -1e-6 * np.cos(2 * np.pi * x)

def vz_(x, y):
    return 0.0

def Bx_(x, y):
    return 1.0

def By_(x, y):
    return 1e-6 * np.cos(2 * np.pi * x)

def Bz_(x, y):
    return 0.0

def P_(x, y):
    return 0.1

x = np.arange(nx) * Dx
y = np.arange(ny) * Dy

rho = np.full((ny, nx), rho_(x, y))
vx = np.full((ny, nx), vx_(x, y))
vy = np.full((ny, nx), vy_(x, y))
vz = np.full((ny, nx), vz_(x, y))
Bx = np.full((ny, nx), Bx_(x, y))
By = np.full((ny, nx), By_(x, y))
Bz = np.full((ny, nx), Bz_(x, y))
P = np.full((ny, nx), P_(x, y))