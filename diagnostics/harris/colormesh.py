import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def read_data(file_path):
    data = np.loadtxt(file_path)
    return data

def reshape_data(data, ny, nx):
    reshaped_data = data.reshape((ny, nx, -1), order='F')
    return reshaped_data

nx = 250
ny = 250

file_path = "harrisres/PRK2_15_1801638342.txt"

data = read_data(file_path)

reshaped_data = reshape_data(data, nx, ny)

dx = 0.1
dy = 0.1

# Extract Bx and By
rho = reshaped_data[:, :, 0]
Bx = reshaped_data[:, :, 4]
By = reshaped_data[:, :, 5]
Bz = reshaped_data[:, :, 6]

# Calculate the derivatives
dBy_dx = np.gradient(By, dx, axis=0)
dBx_dy = np.gradient(Bx, dy, axis=1)

# Calculate Jz
Jz = (dBy_dx - dBx_dy)

toPlot = rho

data_min = np.min(toPlot[-1])
data_max = np.max(toPlot[-1]) + 0.2

Norm = Normalize(vmin=data_min, vmax=data_max)

plt.pcolormesh(toPlot.T, cmap='coolwarm', norm = Norm)  
plt.colorbar()
plt.title('Contour Plot of rho at t=15')
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig('harris_rho_t_15.png')
plt.show()
