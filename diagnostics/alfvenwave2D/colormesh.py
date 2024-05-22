import numpy as np
import matplotlib.pyplot as plt
import os

# Function to read data from file
def read_data(file_path):
    data = np.loadtxt(file_path)
    return data

# Function to reshape data
def reshape_data(data, nx, ny):
    reshaped_data = data.reshape((ny, nx, -1), order='F') # Assuming Fortran-like order
    return reshaped_data

# Function to read times
def read_times(file_paths):
    times = []
    for file_path in file_paths:
        time_str = file_path.split('_')[1].split('.')[0].replace('_', '.')
        time = float(time_str)
        times.append(time)
    return times

results_dir = "results/"
file_paths = [results_dir + file for file in os.listdir(results_dir) if file.startswith("URK2_") and file.endswith(".txt")]

nx = 200
ny = 100

data = [read_data(file_path) for file_path in file_paths]
times = read_times(file_paths)

reshaped_data = [reshape_data(d, nx, ny) for d in data]

time_index = 0

qty = reshaped_data[time_index][:, :, 5]

t=times[time_index]

plt.pcolormesh(qty.T, cmap='coolwarm')  
plt.title('Contour Plot of qty')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()