import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import shutil

nx = 128
Dx = 0.1


quantity_name = 'By'
fixed_index = 0
ny = 1

column_names = ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz', 'P']

def read_file(filename):
    df = pd.read_csv(filename, delim_whitespace=True, header=None, names=column_names)
    return df

results_dir = 'whislerwaveres/'
results_dir2 = 'whislerwaveHLL/'

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

"""quantitieshll = {
    'rho': [],
    'vx': [],
    'vy': [],
    'vz': [],
    'Bx': [],
    'By': [],
    'Bz': [],
    'P': []
}
timeshll = []"""

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

"""for filename in os.listdir(results_dir2):
    if filename.startswith("PRK2_") and filename.endswith(".txt"):
        time_str = filename.split('_')[1]+'.'+filename.split('_')[2].split('.')[0]
        time = float(time_str)
        timeshll.append(time)
        
        df = read_file(os.path.join(results_dir2, filename))

        for quantity in quantitieshll.keys():
            quantitieshll[quantity].append(df[quantity].values.reshape((ny, nx)))

for quantity in quantitieshll.keys():
    quantitieshll[quantity] = np.array(quantitieshll[quantity])"""

x=Dx*np.arange(nx) + 0.5*Dx

output_dir = 'frames'
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.makedirs(output_dir, exist_ok=True)

def update(frame):
    lx = nx*Dx
    m = 4#int(nx/4)

    k = 2 * np.pi / lx * m

    w = (k**2 /2) *(np.sqrt(1+4/k**2) + 1)
    
    expected_value = 1e-7 * np.cos(k * x -w * times[frame] + 0.5488135)

    plt.clf()
    plt.plot(x, quantities[quantity_name][frame, fixed_index, :], color='blue', marker = 'x', markersize=3) # t,y,x
    #plt.plot(x, quantitieshll[quantity_name][frame, fixed_index, :], color='green', marker = 'o', markersize=3) # t,y,x
    plt.plot(x, expected_value)
    plt.title(f'{quantity_name} at t={times[frame]}')  # Format time to one decimal place
    plt.xlabel('x')
    plt.ylabel(quantity_name)
    plt.grid(True)
    #plt.yscale("log")
    plt.tight_layout()
    
    eps = 1e-7
    min_val = np.min(quantities[quantity_name][:, fixed_index, :]) - eps
    max_val = np.max(quantities[quantity_name][:, fixed_index, :]) + eps
    
    plt.ylim(min_val, max_val)
    plt.savefig(f'{output_dir}/frame_{frame:04d}.png')
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
