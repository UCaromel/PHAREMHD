import numpy as np
import matplotlib.pyplot as plt
import os

nx=100
ny=100

def read_file(filename):
    data = np.loadtxt(filename)
    return data

results_dir = 'src/build/results'

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
        
        data = read_file(os.path.join(results_dir, filename))
        
        for i, quantity in enumerate(quantities.keys()):
            quantities[quantity].append(data[:, i].reshape((ny, nx)))

for quantity in quantities.keys():
    quantities[quantity] = np.array(quantities[quantity])
times = np.array(times)

for quantity in quantities.keys():
    quantities[quantity] = np.stack(quantities[quantity], axis=-1)


rho_vs_x = quantities['rho'][:, :, 0][0, :]

plt.plot(rho_vs_x)
plt.xlabel('x')
plt.ylabel('rho')
plt.title('rho vs x at time index {}'.format(0))
plt.show()