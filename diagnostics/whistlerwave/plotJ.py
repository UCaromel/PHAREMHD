import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

nx = 128
Dx = 0.05

Dt = 0.000625


Jy = []
Jz = []
times = []
results_dir = 'whislerwaveres/'

def read_file(filename):
    df = pd.read_csv(filename, delim_whitespace=True, header=None)
    return df

for filename in os.listdir(results_dir):
    if filename.startswith("Jy_") and filename.endswith(".txt"):
        time_str = filename.split('_')[1].split('.')[0]
        time = float(time_str)
        times.append(time)

        df = read_file(os.path.join(results_dir, filename))
        Jy.append(df.values.reshape((7,135)))

    if filename.startswith("Jz_") and filename.endswith(".txt"):
        
        df = read_file(os.path.join(results_dir, filename))
        Jz.append(df.values.reshape((8,135)))

x=Dx*np.arange(135)-2*Dx

def update(frame):

    m = 10
    lx = nx*Dx
    k = 2*np.pi/lx * m
    expectedJy = -k * np.cos(k*x + k**2 * times[frame]*Dt + 0.5488135)*0.01
    expectedJz = -k * np.sin(k*x + k**2 * times[frame]*Dt + 0.5488135)*0.01

    plt.clf()

    #plt.plot(x, expectedJy, 'y-', marker = '+')
    #plt.plot(x, expectedJz, 'g-', marker = 'x')

    plt.plot(x, Jy[frame][3,:], 'b-', marker = 'o')
    #plt.plot(x, Jz[frame][2,:], 'r--',marker = '*')

    #plt.plot(x, Jy2[frame][2,:], 'y-', marker = 'o')
    #plt.plot(x, Jz2[frame][2,:], 'g--',marker = '*')

    plt.xlabel('x')
    plt.grid(True)
    #plt.yscale("log")
    plt.tight_layout()
    
    eps = 1e-7
    min_val = np.min(Jy) - eps
    max_val = np.max(Jz) + eps
    
    plt.ylim(min_val, max_val)
    plt.axvline(x[2], ls='--', color='k')
    plt.axvline(x[-3], ls='--', color='k')


fig = plt.figure(figsize=(8, 6))
ani = FuncAnimation(fig, update, frames=len(times), interval=500)

# Define a function for manually stepping through the frames
def onclick(event):
    if event.key == 'right':
        ani.frame_seq = ani.new_frame_seq()
        fig.canvas.draw()

# Connect the key press event to the function
fig.canvas.mpl_connect('key_press_event', onclick)

plt.show()

"""plt.plot(x, Jy[10][2,:], 'b-', marker = 'o')
plt.plot(x, Jz[10][2,:], 'r--',marker = '*')
plt.axvline(x[2], ls='--', color='k')
plt.axvline(x[-3], ls='--', color='k')

plt.show()"""