import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

nx = 100
ny = 100

column_names = ['rho', 'rhovx', 'rhovy', 'rhovz', 'Bx', 'By', 'Bz', 'Etot']

def read_file(filename):
    df = pd.read_csv(filename, delim_whitespace=True, header=None, names=column_names)
    return df

results_dir = 'src/results'

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

for filename in os.listdir(results_dir):
    if filename.startswith("URK2_") and filename.endswith(".txt"):
        time = int(filename.split('_')[1].split('.')[0])
        times.append(time)
        
        df = read_file(os.path.join(results_dir, filename))

        for quantity in quantities.keys():
            quantities[quantity].append(df[quantity].values.reshape((ny, nx)))

for quantity in quantities.keys():
    quantities[quantity] = np.array(quantities[quantity])
times = np.array(times)


quantity_name = 'By'
fixed_index = 0

def update(frame):
    plt.clf()
    plt.plot(quantities[quantity_name][frame, fixed_index, :], color='blue') # t,y,x
    plt.title(f'{quantity_name} at y={fixed_index}, t={times[frame]}')
    plt.xlabel('x')
    plt.ylabel(quantity_name)
    plt.grid(True)
    plt.tight_layout()
    
    min_val = np.min(quantities[quantity_name][:, fixed_index, :])
    max_val = np.max(quantities[quantity_name][:, fixed_index, :])
    
    plt.ylim(min_val, max_val)


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
