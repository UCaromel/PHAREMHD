import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

nx = 1000
ny = 10
Dt=0.008

quantity_name = 'By'
fixed_index = 0

lx=10

column_names = ['rho', 'rhovx', 'rhovy', 'rhovz', 'Bx', 'By', 'Bz', 'Etot']

def read_file(filename):
    df = pd.read_csv(filename, delim_whitespace=True, header=None, names=column_names)
    return df

results_dir = './results'

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
        # Extract time from filename and format it properly
        time_str = filename.split('_')[1].split('.')[0]
        time_str = time_str.replace('_', '.')  # Replace underscore with dot
        time = float(time_str)
        times.append(time)
        
        df = read_file(os.path.join(results_dir, filename))

        for quantity in quantities.keys():
            quantities[quantity].append(df[quantity].values.reshape((ny, nx)))

for quantity in quantities.keys():
    quantities[quantity] = np.array(quantities[quantity])


times = np.array(times)
times=times*Dt


Dx=lx/nx
x=Dx*np.arange(nx)


def update(frame):    
    t=times[frame]
    expected=1e-6*np.cos(2*np.pi*(x-t)/10)

    plt.clf()
    plt.plot(x,quantities[quantity_name][frame, fixed_index, :], color='blue') # t,y,x
    plt.plot(x,quantities["rhovx"][frame, fixed_index, :], color='red') # t,y,x
    plt.plot(x,expected)
    plt.title(f'{quantity_name} at y={fixed_index}, t={frame}')  # Format time to one decimal place
    plt.xlabel('x')
    plt.ylabel(quantity_name)
    plt.grid(True)
    #plt.yscale("log")
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

def compare_quantities(quantity_name, time1, time2):
    try:
        data1 = quantities[quantity_name][time1, :, :]
        data2 = quantities[quantity_name][time2, :, :]
        if np.array_equal(data1, data2):
            print(f"The data for {quantity_name} at time {time1} is equal to the data at time {time2}")
        else:
            print(f"The data for {quantity_name} at time {time1} is not equal to the data at time {time2}")
    except KeyError:
        print("Invalid quantity name or time")

#compare_quantities('By', 0, 100)