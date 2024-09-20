import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import shutil

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

# for filename in os.listdir(results_dir2):
#     if filename.startswith("PRK2_") and filename.endswith(".txt"):
#         time_str = filename.split('_')[1]+'.'+filename.split('_')[2].split('.')[0]
#         time = float(time_str)
#         times2.append(time)
        
#         df = read_file(os.path.join(results_dir2, filename))

#         for quantity in quantities2.keys():
#             quantities2[quantity].append(df[quantity].values.reshape((ny, nx)))

# for quantity in quantities2.keys():
#     quantities2[quantity] = np.array(quantities2[quantity])

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
    
    expected_value = -1e-2 * np.sin(k * x - w * times[frame] + 0.5488135)

    plt.clf()
    plt.plot(x, quantities[quantity_name][frame, fixed_index, :], color='blue', marker = 'x', markersize=3) # t,y,x
    #plt.plot(x, quantities2[quantity_name][frame, fixed_index, :], color='green', marker = 'o', markersize=3) # t,y,x
    plt.plot(x, expected_value)
    plt.title(f'{quantity_name} at t={times[frame]}')  # Format time to one decimal place
    plt.xlabel('x')
    plt.ylabel(quantity_name)
    plt.grid(True)
    #plt.yscale("log")
    plt.tight_layout()
    
    eps = 1e-3
    min_val = np.min(quantities[quantity_name][:, fixed_index, :]) - eps
    max_val = np.max(quantities[quantity_name][:, fixed_index, :]) + eps
    
    plt.ylim(min_val, max_val)
    plt.savefig(f'{output_dir}/frame_{frame:04d}.png')
    #plt.axvline(x[2]-Dx/2, ls='--', color='k')
    #plt.axvline(x[-2]-Dx/2, ls='--', color='k')

fig = plt.figure(figsize=(8, 6))
ani = FuncAnimation(fig, update, frames=len(times), interval=10)

# Define a function for manually stepping through the frames
def onclick(event):
    if event.key == 'right':
        ani.frame_seq = ani.new_frame_seq()
        fig.canvas.draw()

# Connect the key press event to the function
fig.canvas.mpl_connect('key_press_event', onclick)

plt.show()


# lx = nx*Dx
# m = 4#int(nx/4)

# k = 2 * np.pi / lx * m
# plt.plot(x, quantities[quantity_name][0, fixed_index, :], color='blue', marker = 'x', markersize=3) # t,y,x
# plt.plot(x, 1e-3 * np.cos(k * x + k**2 * times[0] + 0.5488135))
# plt.axvline(x[1]-Dx/2, ls='--', color='k')
# plt.axvline(x[-1]-Dx/2, ls='--', color='k')
# plt.show()

# def calculate_total_energy(quantities, t, gamma):
#     rho = quantities['rho'][t, fixed_index, :]
#     vx = quantities['vx'][t, fixed_index, :]
#     vy = quantities['vy'][t, fixed_index, :]
#     vz = quantities['vz'][t, fixed_index, :]
#     Bx = quantities['Bx'][t, fixed_index, :]
#     By = quantities['By'][t, fixed_index, :]
#     Bz = quantities['Bz'][t, fixed_index, :]
#     P = quantities['P'][t, fixed_index, :]

#     # Calculate kinetic energy density
#     kinetic_energy = 0.5 * rho * (vx**2 + vy**2 + vz**2)
    
#     # Calculate magnetic energy density
#     magnetic_energy = 0.5 * (Bx**2 + By**2 + Bz**2)
    
#     # Calculate internal energy density using pressure and gamma
#     internal_energy = P / (gamma - 1)
    
#     # Total energy density
#     total_energy = internal_energy + kinetic_energy + magnetic_energy
    
#     return np.sum(total_energy)

# print(calculate_total_energy(quantities, 0, 5.0/3.0))
# print(calculate_total_energy(quantities, -1, 5.0/3.0))
