import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

nx = 800
ny = 1

quantities_to_plot = ['rho', 'vx', 'vy', 'By', 'P']
fixed_index = 0
lx = 1

column_names = ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'P']

def read_file(filename):
    df = pd.read_csv(filename, delim_whitespace=True, header=None, names=column_names)
    return df

results_dir = 'shockres'

quantities = {quantity: [] for quantity in column_names}
times = []

for filename in os.listdir(results_dir):
    if filename.startswith("PRK2_") and filename.endswith(".txt"):
        # Extract time from filename and format it properly
        time_str = filename.split('_')[2].split('.')[0]
        time_str = time_str.replace('_', '.')  # Replace underscore with dot
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

# Determine the range for x-axis
x_range = (x.min(), x.max())

# Create figure and grid specification
fig = plt.figure(figsize=(14, 14))
gs = plt.GridSpec(3, 2)

# Plot 'rho' on the centered subplot spanning both columns
ax_rho = plt.subplot(gs[0, :])
ax_rho.plot(x, quantities['rho'][-1, fixed_index, :], 'o', color='blue', markersize=3)
ax_rho.set_title(f'rho at y={fixed_index}, t={times[-1]:.2f}')
ax_rho.set_xlabel('x')
ax_rho.set_ylabel('rho')
ax_rho.grid(True)
eps = 0.1
min_val = np.min(quantities['rho'][:, fixed_index, :]) - eps
max_val = np.max(quantities['rho'][:, fixed_index, :]) + eps
ax_rho.set_xlim(x_range)
ax_rho.set_ylim(min_val, max_val)

# Ensure square aspect ratio for 'rho' subplot
ax_rho.set_aspect((x_range[1]-x_range[0])/(max_val-min_val))

# Plot the remaining quantities in 2x2 configuration
for i, quantity in enumerate(quantities_to_plot[1:]):
    row = (i // 2) + 1
    col = i % 2
    ax = plt.subplot(gs[row, col])
    ax.plot(x, quantities[quantity][-1, fixed_index, :], 'o', color='blue', markersize=3)
    ax.set_title(f'{quantity} at y={fixed_index}, t={times[-1]:.2f}')
    ax.set_xlabel('x')
    ax.set_ylabel(quantity)
    ax.grid(True)
    min_val = np.min(quantities[quantity][:, fixed_index, :]) - eps
    max_val = np.max(quantities[quantity][:, fixed_index, :]) + eps
    ax.set_xlim(x_range)
    ax.set_ylim(min_val, max_val)
    
    # Ensure square aspect ratio for each subplot
    ax.set_aspect((x_range[1]-x_range[0])/(max_val-min_val))

plt.tight_layout()
plt.show()
