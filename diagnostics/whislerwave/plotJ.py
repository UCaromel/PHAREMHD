import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os


Dx = 0.2

Jz = []
Jy = []
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
        Jy.append(df.values.reshape((5,37)))

    if filename.startswith("Jz_") and filename.endswith(".txt"):
        
        df = read_file(os.path.join(results_dir, filename))
        Jz.append(df.values.reshape((6,37)))

x=Dx*np.arange(37)-2*Dx

def update(frame):
    plt.clf()
    plt.plot(x, Jy[frame][2,:], 'b-', marker = 'o')
    plt.plot(x, Jz[frame][2,:], 'r--',marker = '*')
    plt.xlabel('x')
    plt.grid(True)
    #plt.yscale("log")
    plt.tight_layout()
    
    eps = 0.01
    min_val = np.min(Jy) - eps
    max_val = np.max(Jz) + eps
    
    plt.ylim(min_val, max_val)


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

plt.plot(x, Jy[10][2,:], 'b-', marker = 'o')
plt.plot(x, Jz[10][2,:], 'r--',marker = '*')
plt.axvline(x[2], ls='--', color='k')
plt.axvline(x[-3], ls='--', color='k')

plt.show()