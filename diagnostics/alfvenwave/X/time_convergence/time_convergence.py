import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

nx = 10000
ny = 3
Dx = 0.0001

#time_index = 0
studied_quantity = 'By'

results_dir = './time_results'

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

Dt_values = []  
times = {}  # Changed to store times for each stepDt
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
    
    # Initialize the time dictionary entry for this stepDt if it doesn't exist
    if stepDt not in times:
        times[stepDt] = []
    
    # Extract time from filename and format it properly
    time_str = filename.split('_')[2].split('.')[0]  # Extract the part after time integrator name
    time_str = time_str.replace('_', '.')  # Replace underscore with dot
    time = float(time_str)
    times[stepDt].append(time)

    df = pd.read_csv(os.path.join(results_dir, filename), delim_whitespace=True, header=None, names=['rho', 'rhovx', 'rhovy', 'rhovz', 'Bx', 'By', 'Bz', 'Etot'])
    for quantity in quantities[integrator_name]:
        quantities[integrator_name][quantity][stepDt].append(df[quantity].values.reshape((ny, nx)))   # Store the entire data array for each quantity

# Convert lists of times to numpy arrays and remove duplicates
for integrator_name in time_integrators:
    for stepDt in quantities[integrator_name][studied_quantity]:
        times[stepDt] = np.array(times[stepDt])
        times[stepDt] = np.unique(times[stepDt])

        factor = 2 ** stepDt
        
        x = np.arange(nx)*Dx
        expected_value = 1e-6 * np.cos(2 * np.pi * (x - times[stepDt][-1] * 0.00005 / factor))
        
        # Extract computed value for this time
        computed_value = quantities[integrator_name][studied_quantity][stepDt][-1][0, :]
        
        # Calculate error for this Dt
        error = calculate_error(computed_value, expected_value)

        errors[integrator_name].append(error)

print(times[0])

Dt_values = [0.00005 / (2 ** i) for i in range(len(quantities[time_integrators[0]][studied_quantity]))]

plt.figure(figsize=(12, 8))
for integrator_name in time_integrators:
    plt.loglog(Dt_values, errors[integrator_name], marker='o', label=integrator_name)
plt.xlabel('Dt')
plt.ylabel('L2 Norm of Error')
plt.title('L2 Norm of Error vs. Dt')
plt.legend()
plt.grid(True)
plt.show()

expected_value = 1e-6 * np.cos(2 * np.pi * (x - times[0][-1] * 0.00005))
for integrator_name in time_integrators:
    for stepDt in quantities[integrator_name][studied_quantity]:
            
        x = np.arange(nx)*Dx
            
        computed_value = quantities[integrator_name][studied_quantity][stepDt][-1][0, :]


        plt.plot(x, computed_value, label = f"dt = {Dt_values[stepDt]}")
    plt.plot(x, expected_value, label = 'expected')
    plt.legend()
    plt.show()
