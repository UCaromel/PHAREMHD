import numpy as np
import matplotlib.pyplot as plt

def read_data(file_path):
    data = np.loadtxt(file_path)
    return data

def reshape_data(data, ny, nx):
    reshaped_data = data.reshape((ny, nx, -1), order='F')
    return reshaped_data

nx = 100
ny = 100

file_path = "MHDrotorres/URK2_0_0.txt"

data = read_data(file_path)

reshaped_data = reshape_data(data, nx, ny)

rho = reshaped_data[:, :, 0]

plt.pcolormesh(rho.T, cmap='coolwarm')  
plt.title('Contour Plot of rho at t=0')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()