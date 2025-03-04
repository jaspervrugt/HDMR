# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE    222222  11 #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE        22 22   11 #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       22    11 #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE         22     11 #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE    222222  11 #
#                                                                         #
# ####################################################################### #

import sys
import os
import numpy as np

# Get the current working directory
current_directory = os.getcwd()
# Go up one directory
parent_directory = os.path.abspath(os.path.join(current_directory, '..'))
# add this to path
sys.path.append(parent_directory)
# Add another directory
misc_directory = os.path.abspath(os.path.join(parent_directory, 'miscellaneous'))
# add this to path
sys.path.append(misc_directory)

from HDMR import HDMR
from evaluate_TE import evaluate_TE

d = 4                   # Number of parameters
N = 500                 # Number of samples to be used
solver = 'explicit'     # Numerical solution: 'explicit' or 'implicit'

# Check if the file exists
filename = f'XY_{solver}.mat'
try:
    # Try loading the existing data
    data = sio.loadmat(filename)
    X = data['X']
    Y = data['Y']
    print(f"Loaded data from {filename}")
except FileNotFoundError:
    # If file does not exist, generate the data and save it
    X, label_par, Y = evaluate_TE(N, solver)
    sio.savemat(filename, {'X': X, 'Y': Y})
    print(f"Data generated and saved as {filename}")

# Set up options for HDMR toolbox
options = {
    'graphics': 0,
    'maxorder': 2,
    'maxiter': 100,
    'bf1': 1,
    'bf2': 0,
    'bf3': 0,
    'm': 5,
    'K': 50,
    'R': 200,
    'method': 1,
    'alfa': 0.01,
    'lambda': 0.10,
    'vartol': 1e-3,
    'refit': 1
}

# Number of simulations and values per trial
N, n = Y.shape

# Initialize results
results = np.empty((5, 13, n), dtype=object)

if __name__ == '__main__':
	# Run HDMR for each simulation time (temperature)
	for t in range(n):
    		S, Ss, Fx, Em, Xy, RT = HDMR(X, Y[:, t], options)
		results[0:5, 0:13, t] = [Ss, Fx, Em, Xy]

	# Extract total sensitivity of x1, x2, x1x2, and sum
	out = results[1:5, 6, :].reshape(4, 201)

	# Plot sensitivities as a function of temperature
	T = np.arange(500, 701)  # Temperatures from 500 to 700

	plt.plot(T, out[0, :], 'r', linewidth=2, label='$x_{1}$')
	plt.plot(T, out[1, :], 'b', linewidth=2, label='$x_{2}$')
	plt.plot(T, out[2, :], 'g', linewidth=2, label='$x_{1}/x_{2}$')
	plt.plot(T, out[3, :], 'k', linewidth=2, label='$\sum$')

	# Set up the legend and labels
	plt.legend(fontsize=18, loc='best', frameon=True)
	plt.xlabel('Temperature', fontsize=18, fontweight='bold')
	plt.ylabel('Sensitivity', fontsize=18, fontweight='bold')
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.show()