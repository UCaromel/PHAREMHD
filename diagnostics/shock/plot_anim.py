import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import shutil

nx = 800
ny = 1

quantity_name = 'rho'
fixed_index = 0

lx = 1

column_names = ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'P']

def read_file(filename):
    df = pd.read_csv(filename, delim_whitespace=True, header=None, names=column_names)
    return df

results_dir = 'shockres'

quantities = {
    'rho': [],
    'vx': [],
    'vy': [],
    'vz': [],
    'Bx': [],
    'By': [],
    'Bz': [],
    'P': []
}
times = []

# Load data
for filename in os.listdir(results_dir):
    if filename.startswith("PRK2_") and filename.endswith(".txt"):
        time_str = filename.split('_')[1]+'.'+filename.split('_')[2].split('.')[0]
        time = float(time_str)
        times.append(time)
        
        df = read_file(os.path.join(results_dir, filename))

        for quantity in quantities.keys():
            quantities[quantity].append(df[quantity].values.reshape((ny, nx)))

for quantity in quantities.keys():
    quantities[quantity] = np.array(quantities[quantity])

times = np.array(times)

Dx = lx / nx
x = Dx * np.arange(nx) + 0.5 * Dx

# Set up directory for saving frames
output_dir = 'frames'
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.makedirs(output_dir, exist_ok=True)

# Create figure
fig, ax = plt.subplots(figsize=(8, 6))

def update(frame):
    ax.clear()  # Clear the previous frame
    ax.plot(x, quantities[quantity_name][frame, fixed_index, :], 'o', color='blue', markersize=3)
    ax.set_title(f'{quantity_name} at y={fixed_index}, t={times[frame]:.2f}')
    ax.set_xlabel('x')
    ax.set_ylabel(quantity_name)
    ax.grid(True)
    
    eps = 0.1
    min_val = np.min(quantities[quantity_name][:, fixed_index, :]) - eps
    max_val = np.max(quantities[quantity_name][:, fixed_index, :]) + eps
    ax.set_ylim(min_val, max_val)

    # Save the current frame as an image
    plt.tight_layout()
    plt.savefig(f'{output_dir}/frame_{frame:04d}.png')

ani = FuncAnimation(fig, update, frames=len(times), interval=100)

# Show the plot (optional if you just want to save the frames)
plt.show()
