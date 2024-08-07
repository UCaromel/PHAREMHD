import numpy as np
import matplotlib.pyplot as plt

def read_data(file_path):
    data = np.loadtxt(file_path)
    return data

def reshape_data(data, ny, nx):
    reshaped_data = data.reshape((ny, nx, -1), order='F')
    return reshaped_data

nx = 250
ny = 250

file_path = "hallharrisresCT/PRK2_0_5002794389.txt"

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

toPlot = Bz

plt.pcolormesh(toPlot.T, cmap='coolwarm')  
plt.title('Contour Plot of Bz at t=0.5')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

# x_middle = nx // 2

# Jz_cutx = Jz[x_middle, :]

# y_positions = np.arange(ny)

# # Create the plot
# plt.figure(figsize=(10, 6))
# plt.plot(y_positions, Jz_cutx, marker='o')
# plt.title(f'1D cut of Jz at X = {x_middle/nx}')
# plt.xlabel('y')
# plt.ylabel('Jz')
# plt.grid(True)

# plt.show()
