import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

results_dir = './space_results'

nx_initial = 32
ny = 1

fixed_index = 0
time_index = -1

studied_quantity = 'By'

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
n_errors = 5
times = []
errors = np.zeros(n_errors)
nx_values = {}

# List and sort the files in the results directory
filenames = os.listdir(results_dir)
filenames.sort(key=lambda x: (int(x.split('_')[0]), float(x.split('_')[2].split('.')[0].replace('_', '.'))))

for filename in filenames:
    stepDx_str = filename.split('_')[0]  # Extract the part before "URK2_"
    stepDx = int(stepDx_str)
    
    # Compute nx for this stepDx
    nx = nx_initial * (2 ** stepDx)
    nx_values[stepDx] = nx
    
    # Create a new dictionary entry for this Dx if it doesn't exist
    if stepDx not in quantities['rho']:
        for key in quantities:
            quantities[key][stepDx] = []
    
    # Extract time from filename and format it properly
    time_str = filename.split('_')[2]+'.'+filename.split('_')[3].split('.')[0] # Extract the part after "URK2_" (length of "URK2_" is 5)
    time_str = time_str.replace('_', '.')  # Replace underscore with dot
    time = float(time_str)
    times.append(time)

    df = pd.read_csv(os.path.join(results_dir, filename), delim_whitespace=True, header=None, names=['rho', 'rhovx', 'rhovy', 'rhovz', 'Bx', 'By', 'Bz', 'Etot'])
    for quantity in quantities:
        quantities[quantity][stepDx].append(df[quantity].values.reshape(nx))   # Store the entire data array for each quantity

times = np.array(times)
times = np.unique(times)

Dx =np.asarray( [0.4 / (2 ** i) for i in range(len(nx_values))])

# Function to calculate L2 norm of error
def calculate_error(computed, expected, nx):
    #return np.linalg.norm(computed -expected)
    return np.sum(abs(computed - expected))/nx

for stepDx in quantities[studied_quantity]:
    # Get nx for this stepDx
    nx = nx_values[stepDx]
    
    x = Dx[stepDx] * np.arange(nx) + 0.5*Dx[stepDx]

    lx = nx*Dx[stepDx]
    m = int(nx/4)

    k = 2 * np.pi / lx * m
    w = (k**2 /2) *(np.sqrt(1+4/k**2) + 1)
    
    expected_value = 1e-7 * np.cos(k * x - w * times[time_index] + 0.5488135)
    
    # Extract computed value for this time_index
    computed_value = quantities[studied_quantity][stepDx][time_index]

    #print(computed_value, expected_value, "#####################")
    
    # Calculate error for this Dx
    error = calculate_error(computed_value, expected_value , nx) #/ error0
    
    errors[stepDx] = error

# Assuming Dx and errors are your data arrays
# Log-transform the data
log_Dx = np.log(Dx)
log_errors = np.log(errors)

# Perform a linear fit to the log-log data
coefficients = np.polyfit(log_Dx, log_errors, 1)

# Extract the slope (order of accuracy) and intercept
slope, intercept = coefficients

# Generate the fitted line for plotting
fitted_log_errors = np.polyval(coefficients, log_Dx)
fitted_errors = np.exp(fitted_log_errors)

# Plot the original data and the fitted line
plt.figure(figsize=(12, 8))
plt.loglog(Dx, errors, marker='o', label='Numerical Error')
#plt.loglog(Dx,np.exp(log_Dx*slope + intercept), label="manual")
plt.loglog(Dx, fitted_errors, linestyle='--', label=f'Fitted Line (slope={slope:.2f})')
plt.xlabel('Dx')
plt.ylabel('L2 Norm of Error')
plt.title('L2 Norm of Error vs. Dx')
plt.grid(True)
plt.legend()
plt.show()

for stepDx in quantities[studied_quantity]:
    computed_value = quantities[studied_quantity][stepDx][time_index]
    plt.plot(Dx[stepDx]*np.arange(nx_values[stepDx]), computed_value, label=f"dx = {Dx[stepDx]}")
plt.legend()    
