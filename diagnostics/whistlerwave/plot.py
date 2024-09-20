import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

nx = 128
Dx = 0.05


quantity_name = 'Bz'
fixed_index = 0
ny = 1

column_names = ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'P']

def read_file(filename):
    df = pd.read_csv(filename, delim_whitespace=True, header=None, names=column_names)
    return df

results_dir = 'whislerwaveres_singlemode/'

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

# Calculate w
lx = nx * Dx
m = 4
k = 2 * np.pi / lx * m
w = (k**2 / 2) * (np.sqrt(1 + 4/k**2) + 1)

# Calculate the desired time points
t0 = 0
t1 = (2 * np.pi / w) / 2
t2 = 2 * np.pi / w

# Find the closest time indices
time_indices = [
    np.argmin(np.abs(np.array(times) - t)) 
    for t in [t0, t1, t2]
]

x=Dx*np.arange(nx) + 0.5*Dx

# Create a figure with 3 subplots
fig, axs = plt.subplots(3, 1, figsize=(10, 15))

for i, (ax, time_index) in enumerate(zip(axs, time_indices)):
    t = times[time_index]
    
    expected_value = -1e-2 * np.sin(k * x - w * t + 0.5488135)
    
    ax.plot(x, quantities[quantity_name][time_index, fixed_index, :], color='blue', marker='x', markersize=3, label='Simulated')
    ax.plot(x, expected_value, color='red', label='Expected')
    
    ax.set_title(f'{quantity_name} at t={t:.2f}')
    ax.set_xlabel('x')
    ax.set_ylabel(quantity_name)
    ax.grid(True)
    ax.legend()
    
    eps = 1e-3
    min_val = np.min(quantities[quantity_name][:, fixed_index, :]) - eps
    max_val = np.max(quantities[quantity_name][:, fixed_index, :]) + eps
    ax.set_ylim(min_val, max_val)

plt.tight_layout()
plt.savefig('three_time_points_comparison.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.show()