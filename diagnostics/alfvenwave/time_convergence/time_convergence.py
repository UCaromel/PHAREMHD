import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

nx = 100
ny = 100
Dx = 0.01

time_index = 100
studied_quantity = 'By'

# Define the directory where results are stored
results_dir = './results'

# Define quantities as a dictionary of dictionaries for each time integrator
time_integrators = ['UEuler', 'URK2', 'URK3']

quantities = {integrator: {
    'rho': {},
    'rhovx': {},
    'rhovy': {},
    'rhovz': {},
    'Bx': {},
    'By': {},
    'Bz': {},
    'Etot': {},
} for integrator in time_integrators}

Dt_values = []  # Array to store Dt values
times = {integrator: [] for integrator in time_integrators}
errors = {integrator: [] for integrator in time_integrators}

# Function to calculate L2 norm of error
def calculate_error(computed, expected):
    return np.linalg.norm(computed - expected)

# Loop through the files in the results directory
for filename in os.listdir(results_dir):
    stepDt_str = filename.split('_')[0]  # Extract the part before time integrator name
    stepDt = int(stepDt_str)
    integrator_name = filename.split('_')[1]  # Extract the time integrator name
    
    # Create a new dictionary entry for this Dt if it doesn't exist
    if stepDt not in quantities[integrator_name]['rho']:
        for key in quantities[integrator_name]:
            quantities[integrator_name][key][stepDt] = []
    
    # Extract time from filename and format it properly
    time_str = filename.split('_')[2].split('.')[0]  # Extract the part after time integrator name
    time_str = time_str.replace('_', '.')  # Replace underscore with dot
    time = float(time_str)
    times[integrator_name].append(time)

    df = pd.read_csv(os.path.join(results_dir, filename), delim_whitespace=True, header=None, names=['rho', 'rhovx', 'rhovy', 'rhovz', 'Bx', 'By', 'Bz', 'Etot'])
    for quantity in quantities[integrator_name]:
        quantities[integrator_name][quantity][stepDt].append(df[quantity].values.reshape((ny, nx)))   # Store the entire data array for each quantity




# Adjust the expected value calculation for each successive Dt
for integrator_name in time_integrators:
    times[integrator_name] = np.array(times[integrator_name])
    times[integrator_name] = np.unique(times[integrator_name])
    for stepDt in quantities[integrator_name][studied_quantity]:
        factor = 10 ** stepDt
        
        x = np.arange(nx)*Dx
        expected_value = 1e-6 * np.cos(2 * np.pi * (x - times[integrator_name][time_index] * 0.008 / factor))
        
        # Extract computed value for this time
        computed_value = quantities[integrator_name][studied_quantity][stepDt][time_index][0, :]
        
        # Calculate error for this Dt
        error = calculate_error(computed_value, expected_value)
        
        # Append error to list
        errors[integrator_name].append(error)


Dt_values = [0.008 / (10 ** i) for i in range(len(quantities[time_integrators[0]][studied_quantity]))]  # Adjusted for initial Dt = 0.008

# Plot errors against Dt for each time integrator
plt.figure(figsize=(12, 8))
for integrator_name in time_integrators:
    plt.loglog(Dt_values, errors[integrator_name], marker='o', label=integrator_name)
plt.xlabel('Dt')
plt.ylabel('L2 Norm of Error')
plt.title('L2 Norm of Error vs. Dt')
plt.legend()
plt.grid(True)
plt.show()
