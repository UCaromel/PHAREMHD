import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import os

# Function to read data from file
def read_data(file_path):
    data = np.loadtxt(file_path)
    return data

# Function to reshape data
def reshape_data(data, ny, nx):
    reshaped_data = data.reshape((ny, nx, -1), order='F') # Assuming Fortran-like order
    return reshaped_data

# Function to read times
def read_times(file_paths):
    times = []
    for file_path in file_paths:
        time_str = file_path.split('_')[1]+'.'+file_path.split('_')[2].split('.')[0]
        time = float(time_str)
        times.append(time)
    return times

results_dir = "orszagtangCTAverage/"
file_paths = [results_dir + file for file in os.listdir(results_dir) if file.startswith("PRK2_") and file.endswith(".txt")]

nx = 128
ny = 128

studied_index = 7

data = [read_data(file_path) for file_path in file_paths]
times = read_times(file_paths)

reshaped_data = [reshape_data(d, nx, ny) for d in data]

data_min = np.min(reshaped_data[-1][:, :, studied_index])
data_max = np.max(reshaped_data[-1][:, :, studied_index])

Norm = Normalize(vmin=data_min, vmax=data_max)

qty = reshaped_data[-1][:,:,studied_index]

plt.pcolormesh(qty.T, cmap='coolwarm', norm=Norm)  
plt.title('Contour Plot of qty')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()