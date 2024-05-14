import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

nx = 1000
ny = 10
Dt = 0.8e-6

lx = 10

fixed_index = 0
time_index = 1

# Define the directory where results are stored
results_dir = './results'

# Define quantities as a dictionary of dictionaries
quantities = {
    'rho': {},
    'rhovx': {},
    'rhovy': {},
    'rhovz': {},
    'Bx': {},
    'By': {},
    'Bz': {},
    'Etot': {},
}

Dx_values = []  # Array to store Dx values
times = []
errors = []

# Function to calculate L2 norm of error
def calculate_error(computed, expected):
    return np.linalg.norm(computed - expected)

# Loop through the files in the results directory
for filename in os.listdir(results_dir):
    stepDx_str = filename.split('_')[0]  # Extract the part before "URK2_"
    stepDx = int(stepDx_str)
    
    # Create a new dictionary entry for this Dx if it doesn't exist
    if stepDx not in quantities['rho']:
        for key in quantities:
            quantities[key][stepDx] = []
    
    # Extract time from filename and format it properly
    time_str = filename.split('_')[2].split('.')[0]  # Extract the part after "URK2_" (length of "URK2_" is 5)
    time_str = time_str.replace('_', '.')  # Replace underscore with dot
    time = float(time_str)
    times.append(time)

    df = pd.read_csv(os.path.join(results_dir, filename), delim_whitespace=True, header=None, names=['rho', 'rhovx', 'rhovy', 'rhovz', 'Bx', 'By', 'Bz', 'Etot'])
    for quantity in quantities:
        quantities[quantity][stepDx].append(df[quantity].values.reshape((ny, nx)))   # Store the entire data array for each quantity

times = np.array(times)
times = np.unique(times)

for stepDx in quantities['By']:
    factor = 10 ** (stepDx)
    Dx = 0.01 

    x = float(Dx) * np.arange(nx)
    expected_value = 1e-6 * np.cos(2 * np.pi * (x - times[time_index] * (Dt / factor)) / 10)
    
    # Extract computed value for this time
    computed_value = quantities['By'][stepDx][time_index][fixed_index, :]
    
    # Calculate error for this Dx
    error = calculate_error(computed_value, expected_value)
    
    # Append error to list
    errors.append(error)

Dx = [0.01 / (10 ** i) for i in range(3, -1, -1)]

# Plot errors against time
plt.figure(figsize=(12, 8))
plt.loglog(Dx, errors, marker='o')
plt.xlabel('Time')
plt.ylabel('L2 Norm of Error')
plt.title('L2 Norm of Error vs. Time')
plt.grid(True)
plt.show()

print(errors[2])
