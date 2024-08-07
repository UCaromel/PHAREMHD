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

file_path1 = "hallharrisresCT/PRK2_0_5002794389.txt"
file_path2 = "hallharrisresCT2/PRK2_0_5000859584.txt"
# file_path1 = "hallharrisres2/PRK2_2_4174122015.txt"
# file_path2 = "hallharrisresCT2long/PRK2_2_4039801735.txt"

# Read and process data from the first file
data1 = read_data(file_path1)
reshaped_data1 = reshape_data(data1, ny, nx)

# Extract Bz for the first file
Bz1 = reshaped_data1[:, :, 6]

# Read and process data from the second file
data2 = read_data(file_path2)
reshaped_data2 = reshape_data(data2, ny, nx)

# Extract Bz for the second file
Bz2 = reshaped_data2[:, :, 6]

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Plot Bz for the first file
im1 = axes[0].pcolormesh(Bz1.T, cmap='coolwarm')
axes[0].set_title('Contour Plot of Bz at t=0.5 (NoCT)')
axes[0].set_xlabel('X')
axes[0].set_ylabel('Y')
fig.colorbar(im1, ax=axes[0])

# Plot Bz for the second file
im2 = axes[1].pcolormesh(Bz2.T, cmap='coolwarm')
axes[1].set_title('Contour Plot of Bz at t=0.5 (CT)')
axes[1].set_xlabel('X')
axes[1].set_ylabel('Y')
fig.colorbar(im2, ax=axes[1])

plt.tight_layout()
plt.show()
