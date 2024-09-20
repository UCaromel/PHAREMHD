import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import shutil

# Parameters
nx = 50
ny = 3
Dt = 5e-4
quantity_name = 'By'
fixed_index = 0
lx = 1

# Define column names
column_names = ['rho', 'rhovx', 'rhovy', 'rhovz', 'Bx', 'By', 'Bz', 'Etot']

def read_file(filename):
    df = pd.read_csv(filename, delim_whitespace=True, header=None, names=column_names)
    return df

# Results directory
results_dir = './space_results'

# Dictionary to store quantities
quantities = {
    'rho': [],
    'rhovx': [],
    'rhovy': [],
    'rhovz': [],
    'Bx': [],
    'By': [],
    'Bz': [],
    'Etot': []
}

times = []

# Read all files in the results directory
for filename in os.listdir(results_dir):
    if filename.startswith("0_URK2_") and filename.endswith(".txt"):
        # Extract time from filename and format it properly
        time_str = filename.split('_')[2]+'.'+filename.split('_')[3].split('.')[0]
        time = float(time_str)
        times.append(time)
        
        df = read_file(os.path.join(results_dir, filename))

        for quantity in quantities.keys():
            quantities[quantity].append(df[quantity].values.reshape((ny, nx)))

# Convert lists to numpy arrays
for quantity in quantities.keys():
    quantities[quantity] = np.array(quantities[quantity])

times = np.array(times)
times = times

# Define x-axis
Dx = lx / nx
x = Dx * np.arange(nx) + 0.5 * Dx

# Create output directory for frames
output_dir = 'frames'
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.makedirs(output_dir, exist_ok=True)

# Define the update function for the animation
def update(frame):
    plt.clf()  # Clear the previous plot
    t = times[frame]

    # Plot the desired quantity
    plt.plot(x, quantities[quantity_name][frame, fixed_index, :], color='blue')
    plt.title(f'{quantity_name} at y={fixed_index}, t={times[frame]:.4f}')
    plt.xlabel('x')
    plt.ylabel(quantity_name)
    plt.grid(True)
    plt.tight_layout()
    
    # Set y-limits based on data range
    min_val = np.min(quantities[quantity_name][:, fixed_index, :])
    max_val = np.max(quantities[quantity_name][:, fixed_index, :])
    plt.ylim(min_val, max_val)

    # Save the current frame
    plt.savefig(f'{output_dir}/frame_{frame:04d}.png')

# Create the figure
fig = plt.figure(figsize=(8, 6))

# Create animation
ani = FuncAnimation(fig, update, frames=len(times), interval=100)

# Manually step through the frames with arrow keys
def onclick(event):
    if event.key == 'right':
        ani.frame_seq = ani.new_frame_seq()
        fig.canvas.draw()

# Connect key press event
fig.canvas.mpl_connect('key_press_event', onclick)

# Show plot (optional, can comment out if just saving frames)
plt.show()
