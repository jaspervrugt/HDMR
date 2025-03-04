# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      666666   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE          66       #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       666666   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE          66  66   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      666666   #
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
from user_function import *

# Case study
d = 3           		# Number of parameters
N = 5000        		# Number of samples
mu = np.array([0.5, 0.5, 0.5])  # Mean values for the distribution
s = np.array([0.2, 0.2, 0.18])  # Standard deviations
r12 = 0.6       		# Correlation between the first two variables

# Create covariance matrix
C = create_covariance(s, r12)

# Draw N samples from a multivariate normal distribution
X = np.random.multivariate_normal(mu, C, N)

# Normalize X
X_min = X.min(axis=0)
X_max = X.max(axis=0)
X = (X - X_min) / (X_max - X_min)

# Define user function parameters
a = np.array([1, 2])
b = np.array([2, 3])
c = np.array([3, 1, 2])
e = np.array([1, 2, 2, 3])

# Evaluate user-defined function y = f(X)
y = user_function(X, mu, a, b, c, e)

# Add random noise
sigma2 = np.var(y) / 100
noise = np.random.normal(0, np.sqrt(sigma2), N)
y = y + noise

# Specify HDMR options (assuming you have a Python implementation of HDMR)
options = {
    'graphics': 1,
    'maxorder': 3,
    'maxiter': 100,
    'bf1': 1,
    'bf2': 0,
    'bf3': 0,
    'm': 2,
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

	print("y (with random noise added):", y)
