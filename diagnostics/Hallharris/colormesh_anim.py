import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import Normalize
import shutil

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

results_dir = "hallharrisres/"
file_paths = [results_dir + file for file in os.listdir(results_dir) if file.startswith("PRK2_") and file.endswith(".txt")]

nx = 250
ny = 250

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

output_dir = 'frames'
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.makedirs(output_dir, exist_ok=True)

def update(frame):
    im.set_array(reshaped_data[frame][:, :, studied_index].T)
    plt.title(f'rho at t={times[frame]}')
    #plt.savefig(f'{output_dir}/frame_{frame:04d}.png')
    return im,

ani = FuncAnimation(fig, update, frames=len(times), interval=100)

#gif_writer = PillowWriter(fps=10)  # Adjust fps as needed
#ani.save('hallharris.gif', writer=gif_writer)

plt.xlabel('X')
plt.ylabel('Y')
plt.show()

#ffmpeg -r 10 -i frames/frame_%04d.png -vcodec mpeg4 -q:v 5 hallharris.mp4