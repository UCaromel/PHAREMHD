import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

nx = 1
ny = 128
Dx = 1
Dy = 0.1

nghost = 2

nghostJ = nghost + 1

Dt = 0.04

Jx = []
Jy = []
Jz = []
times = []
results_dir = 'whislerwaveres/'

def read_file(filename):
    df = pd.read_csv(filename, delim_whitespace=True, header=None)
    return df

for filename in os.listdir(results_dir):
    if filename.startswith("Jx_") and filename.endswith(".txt"):
        time_str = filename.split('_')[1].split('.')[0]
        time = float(time_str)
        times.append(time)

        df = read_file(os.path.join(results_dir, filename))
        Jx.append(df.values.reshape((ny + 2*nghostJ + 1, nx + 2*nghostJ)))

    if filename.startswith("Jy_") and filename.endswith(".txt"):
        
        df = read_file(os.path.join(results_dir, filename))
        Jy.append(df.values.reshape((ny + 2*nghostJ, nx + 2*nghostJ + 1)))

    if filename.startswith("Jz_") and filename.endswith(".txt"):
        
        df = read_file(os.path.join(results_dir, filename))
        Jz.append(df.values.reshape((ny + 2*nghostJ + 1, nx + 2*nghostJ + 1)))

x=Dx*np.arange(nx + 2*nghostJ )-nghostJ*Dx
y=Dy*np.arange(ny + 2*nghostJ+1)-nghostJ*Dy

# def update(frame):

#     m = 4
#     ly = ny*Dy
#     k = 2*np.pi/ly * m
#     #expectedJy = -k * np.cos(k*x + k**2 * times[frame]*Dt + 0.5488135)*0.01
#     #expectedJz = -k * np.sin(k*x + k**2 * times[frame]*Dt + 0.5488135)*0.01

#     plt.clf()

#     #plt.plot(x, expectedJy, 'y-', marker = '+')
#     #plt.plot(x, expectedJz, 'g-', marker = 'x')

#     plt.plot(y, Jx[frame][:,3], 'b-', marker = 'o')
#     #plt.plot(x, Jz[frame][2,:], 'r--',marker = '*')

#     #plt.plot(x, Jy2[frame][2,:], 'y-', marker = 'o')
#     #plt.plot(x, Jz2[frame][2,:], 'g--',marker = '*')

#     plt.xlabel('x')
#     plt.grid(True)
#     #plt.yscale("log")
#     plt.tight_layout()
#     plt.axvline(y[2], ls='--', color='k')
#     plt.axvline(y[-3], ls='--', color='k')


# fig = plt.figure(figsize=(8, 6))
# ani = FuncAnimation(fig, update, frames=len(times), interval=500)

# # Define a function for manually stepping through the frames
# def onclick(event):
#     if event.key == 'right':
#         ani.frame_seq = ani.new_frame_seq()
#         fig.canvas.draw()

# # Connect the key press event to the function
# fig.canvas.mpl_connect('key_press_event', onclick)

# plt.show()

plt.plot(y, Jx[7][:,3], 'b-', marker = 'o')
plt.axvline(y[3], ls='--', color='k')
plt.axvline(y[-4], ls='--', color='k')

plt.show()