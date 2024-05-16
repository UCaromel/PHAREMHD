import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

results_dir = './space_results'

ny_initial = 50
nx = 3
Dt = 5e-4

fixed_index = 0
time_index = -1

studied_quantity = 'Bx'

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
ny_values = {}

# List and sort the files in the results directory
filenames = os.listdir(results_dir)
filenames.sort(key=lambda x: (int(x.split('_')[0]), float(x.split('_')[2].split('.')[0].replace('_', '.'))))

for filename in filenames:
    stepDy_str = filename.split('_')[0]  # Extract the part before "URK2_"
    stepDy = int(stepDy_str)
    
    # Compute ny for this stepDy
    ny = ny_initial * (2 ** stepDy)
    ny_values[stepDy] = ny
    
    # Create a new dictionary entry for this Dy if it doesn't exist
    if stepDy not in quantities['rho']:
        for key in quantities:
            quantities[key][stepDy] = []
    
    # Extract time from filename and format it properly
    time_str = filename.split('_')[2].split('.')[0]  # Extract the part after "URK2_" (length of "URK2_" is 5)
    time_str = time_str.replace('_', '.')  # Replace underscore with dot
    time = float(time_str)
    times.append(time)

    df = pd.read_csv(os.path.join(results_dir, filename), delim_whitespace=True, header=None, names=['rho', 'rhovx', 'rhovy', 'rhovz', 'Bx', 'By', 'Bz', 'Etot'])
    for quantity in quantities:
        quantities[quantity][stepDy].append(df[quantity].values.reshape((ny, nx)))   # Store the entire data array for each quantity

times = np.array(times)
times = np.unique(times)

Dy = np.asarray([0.02 / (2 ** i) for i in range(len(ny_values))])

# Function to calculate L2 norm of error
def calculate_error(computed, expected, ny):
    return np.sum(abs(computed - expected)) / ny

for stepDy in quantities[studied_quantity]:
    factor = 2 ** stepDy

    # Get ny for this stepDy
    ny = ny_values[stepDy]
    
    y = Dy[stepDy] * np.arange(ny)

    expected_value = -1e-6 * np.cos(2 * np.pi * (y - times[time_index] * Dt))
    
    # Extract computed value for this time_index
    computed_value = quantities[studied_quantity][stepDy][time_index][:, fixed_index]
    
    # Calculate error for this Dy
    error = calculate_error(computed_value, expected_value, ny)
    
    errors[stepDy] = error

# Log-transform the data
log_Dy = np.log(Dy)
log_errors = np.log(errors)

# Perform a linear fit to the log-log data
coefficients = np.polyfit(log_Dy, log_errors, 1)

# Extract the slope (order of accuracy) and intercept
slope, intercept = coefficients

# Generate the fitted line for plotting
fitted_log_errors = np.polyval(coefficients, log_Dy)
fitted_errors = np.exp(fitted_log_errors)

# Plot the original data and the fitted line
plt.figure(figsize=(12, 8))
plt.loglog(Dy, errors, marker='o', label='Numerical Error')
plt.loglog(Dy, fitted_errors, linestyle='--', label=f'Fitted Line (slope={slope:.2f})')
plt.xlabel('Dy')
plt.ylabel('L2 Norm of Error')
plt.title('L2 Norm of Error vs. Dy')
plt.grid(True)
plt.legend()
plt.show()

# Plot computed values for different Dy
for stepDy in quantities[studied_quantity]:
    computed_value = quantities[studied_quantity][stepDy][time_index][:, fixed_index]
    plt.plot(Dy[stepDy] * np.arange(ny_values[stepDy]), computed_value, label=f"dy = {Dy[stepDy]}")
plt.legend()
plt.xlabel('y')
plt.ylabel(studied_quantity)
plt.title('Computed By values for different Dy')
plt.grid(True)
plt.show()
