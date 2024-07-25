import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
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

results_dir = 'whislerwaveres/'

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
    if filename.startswith("FluxDifx_") and filename.endswith(".txt"):
        time_str = filename.split('_')[1].split('.')[0]
        time = float(time_str)
        times.append(time)
        
        df = read_file(os.path.join(results_dir, filename))

        for quantity in quantities.keys():
            quantities[quantity].append(df[quantity].values.reshape((ny, nx)))

for quantity in quantities.keys():
    quantities[quantity] = np.array(quantities[quantity])

x=Dx*np.arange(nx) + 0.5*Dx

def update(frame):
    
    plt.clf()
    plt.plot(x, quantities[quantity_name][frame, fixed_index, :], color='blue', marker = 'x', markersize=3) # t,y,x
    plt.title(f'{quantity_name} at t={times[frame]}')  # Format time to one decimal place
    plt.xlabel('x')
    plt.ylabel(quantity_name)
    plt.grid(True)
    #plt.yscale("log")
    plt.tight_layout()
    
    eps = 0.001
    min_val = np.min(quantities[quantity_name][:, fixed_index, :]) - eps
    max_val = np.max(quantities[quantity_name][:, fixed_index, :]) + eps
    
    plt.ylim(min_val, max_val)
    #plt.axvline(x[2]-Dx/2, ls='--', color='k')
    #plt.axvline(x[-2]-Dx/2, ls='--', color='k')

fig = plt.figure(figsize=(8, 6))
ani = FuncAnimation(fig, update, frames=len(times), interval=100)

# Define a function for manually stepping through the frames
def onclick(event):
    if event.key == 'right':
        ani.frame_seq = ani.new_frame_seq()
        fig.canvas.draw()

# Connect the key press event to the function
fig.canvas.mpl_connect('key_press_event', onclick)

plt.show()

"""
lx = nx*Dx
m = 4#int(nx/4)

k = 2 * np.pi / lx * m
plt.plot(x, quantities[quantity_name][0, fixed_index, :], color='blue', marker = 'x', markersize=3) # t,y,x
plt.plot(x, 1e-3 * np.cos(k * x + k**2 * times[0] + 0.5488135))
plt.axvline(x[1]-Dx/2, ls='--', color='k')
plt.axvline(x[-1]-Dx/2, ls='--', color='k')
plt.show()
"""
