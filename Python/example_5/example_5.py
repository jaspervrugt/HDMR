# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      555555   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE          55       #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       555555   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE              55   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      555555   #
#                                                                         #
# ####################################################################### #

# Example from the following paper
#  Li, G., and H. Rabitz (2012), General formulation of HDMR component 
#      functions with independent and correlated variables, J. 
#      Math. Chem., 50, pp. 99-130

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

# Parameters
d = 3                                # Number of dimensions
N = 5000                             # Number of samples
mu = 0.5 * np.ones(d)                # Sample mean, µ

# Covariance matrix, Σ (constructed as in the MATLAB code)
C = np.eye(d)                        # Start with an identity matrix
C[0, 1] = 0.6                        # Set C(1, 2) = 0.6
C[1, 0] = 0.6                        # Set C(2, 1) = 0.6

# Draw N samples from N(µ, Σ) - Multivariate normal distribution
X = np.random.multivariate_normal(mu, C, N)

# Normalize X (similar to MATLAB's repmat and normalization)
X_min = np.min(X, axis=0)
X_max = np.max(X, axis=0)
X = (X - X_min) / (X_max - X_min)

# Define the function y = f(x)
y = np.sin(2 * np.pi * X[:, 0] - np.pi) + 7 * (np.sin(2 * np.pi * X[:, 1] - np.pi))**2 \
    + 0.1 * (2 * np.pi * X[:, 2] - np.pi)**4 * np.sin(2 * np.pi * X[:, 0] - np.pi)

# Variance of random error
sigma2 = np.var(y) / 55

# Add random error (Gaussian noise)
y = y + np.random.normal(0, np.sqrt(sigma2), N)

# HDMR options (Python dictionary equivalent to MATLAB's struct)
options = {
    'graphics': 1,
    'maxorder': 3,
    'maxiter': 100,
    'bf1': 1,
    'bf2': 0,
    'bf3': 0,
    'm': 4,
    'K': 10,
    'R': 300,
    'method': 1,
    'alfa': 0.01,
    'lambda': 0.10,
    'vartol': 1e-3,
    'refit': 1
}

if __name__ == '__main__':
	# Run the HDMR toolbox
	S, Ss, Fx, Em, Xy, RT = HDMR(X, y, options)

	# Print the results (optional)
	print("S:", S)
	print("Ss:", Ss)
	print("Fx:", Fx)
	print("Em:", Em)
	print("Xy:", Xy)
