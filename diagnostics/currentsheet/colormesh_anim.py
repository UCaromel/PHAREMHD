import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import Normalize

# Function to read data from file
def read_data(file_path):
    data = np.loadtxt(file_path)
    return data

# Function to reshape data
def reshape_data(data, ny, nx):
    reshaped_data = data.reshape((ny, nx, -1), order='F') 
    return reshaped_data

# Function to read times
def read_times(file_paths):
    times = []
    for file_path in file_paths:
        time_str = file_path.split('_')[1]+'.'+file_path.split('_')[2].split('.')[0]
        time = float(time_str)
        times.append(time)
    return times

results_dir = "currentres/"
file_paths = [results_dir + file for file in os.listdir(results_dir) if file.startswith("PRK2_") and file.endswith(".txt")]

nx = 200
ny = 200

data = [read_data(file_path) for file_path in file_paths]
times = read_times(file_paths)

reshaped_data = [reshape_data(d, nx, ny) for d in data]


studied_index = 0

data_min = np.min(reshaped_data[-1][:, :, studied_index])
data_max = np.max(reshaped_data[-1][:, :, studied_index])

Norm = Normalize(vmin=data_min, vmax=data_max)

fig, ax = plt.subplots()
im = ax.pcolormesh(reshaped_data[0][:, :, studied_index].T, cmap='coolwarm', norm=Norm) 
ax.set_aspect('equal')
fig.colorbar(im, ax=ax)

def update(frame):
    im.set_array(reshaped_data[frame][:, :, studied_index].T)
    plt.title(f'qty at t={times[frame]}')
    return im,

ani = FuncAnimation(fig, update, frames=len(times), interval=100)

#gif_writer = PillowWriter(fps=10)  # Adjust fps as needed
#ani.save('orszagtangP.gif', writer=gif_writer)

plt.xlabel('X')
plt.ylabel('Y')
plt.show()
