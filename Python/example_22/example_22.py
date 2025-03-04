# ####################################################################### #
#                                                                         #
# EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   222222 222222 #
# EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       22 22  22 22  #
# EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE      22     22   #
# EE       XXXX   AAAAAA  MM   MM  PP      LL      EE        22     22    #
# EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   222222 222222 #
#                                                                         #
# ####################################################################### #

# Example of the following paper
#  Gao, Y., A. Sahin, and J.A. Vrugt (2023), Probabilistic Sensitivity 
#      Analysis With Dependent Variables: Covariance-Based Decomposition 
#      of Hydrologic Models, Water Resources Research, 59 (4), 
#      e2022WR0328346, https://doi.org/10.1029/2022WR032834

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
from HDMR_functions import LH_sampling
from rainfall_runoff import rainfall_runoff

# Dimensionality of the model
d = 7
N = 1000  # Number of samples used by HDMR_EXT

# Parname:      Imax Smax Qsmax  alE alF Kfast Kslow
Par_info = {
    'min': np.array([0.5, 10, 0, 1e-6, -10, 0, 0]),  # minimum values
    'max': np.array([10, 1000, 100, 100, 10, 10, 150])  # maximum values
}
label_par = ['$I_{\rm max}$', '$S_{\rm max}$', '$Q_{\rm s,max}$', 
             '$\alpha_{\rm E}$', '$\alpha_{\rm F}$', '$K_{\rm f}$', '$K_{\rm s}$']


## FIX
# Dummy function for rainfall_runoff (replace with actual model)
def rainfall_runoff(x):
    # Placeholder for the actual model function
    return np.sum(x)  # Example, simply summing the parameters

# Generate Latin Hypercube samples for the parameters
X = LH_sampling(Par_info['min'], Par_info['max'], N)

# Initialize Y matrix (based on rainfall_runoff output)
Y = np.zeros((1, N))  # Assuming one output per run for simplicity

# Run the model once to determine n (the output size)
Y[:, 0] = rainfall_runoff(X[0, :d])  # Initial value

# Parallelize the remaining model evaluations
with ProcessPoolExecutor() as executor:
    results = list(executor.map(lambda x: rainfall_runoff(x), X[1:, :d]))

# Fill the Y matrix with results
Y[:, 1:] = np.array(results).T  # Transpose to match the shape

print(Y)

# Specify HDMR structure options as a dictionary
options = {
    'graphics': 0,
    'maxorder': 1,
    'maxiter': 100,
    'bf1': 1,
    'bf2': 0,
    'bf3': 0,
    'm': 5,
    'K': 10,
    'R': 300,
    'method': 1,
    'alfa': 0.01,
    'lambda': 0.10,
    'vartol': 1e-3,
    'refit': 1
}

# Initialize SA as a 3D numpy array (9x15xN), where N is defined elsewhere
ny = Y.shape[0]  # Get the number of rows
N = Y.shape[1]   # Get the number of columns

results = np.empty((9, 15, ny))  # Initialize the results array (adjust size as needed)
if __name__ == '__main__':
	# Run the HDMR toolbox
	for t in range(ny):
    		results[:,:,t], Ss, Fx, Em, XY = HDMR(X, Y[t, :N], options)
