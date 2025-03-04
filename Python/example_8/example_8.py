# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      888888   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE          88  88   #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       888888   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE          88  88   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      888888   #
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
from corr2cov import corr2cov

d = 3    # Number of parameters
N = 5000  # Number of samples

# Mean vector (Sample means µ)
mu = 0.5 * np.ones(d)  # Creates a vector of 0.5s of length d

# Correlation matrix
R = np.array([[1, 0.5, 0.4], 
              [0.5, 1, 0.7], 
              [0.4, 0.7, 1]])

# Sample standard deviations (assuming 1 for each)
std_devs = np.ones(d)

# Covariance matrix
C = corr2cov(std_devs, R)

# Draw N samples from N(mu, C)
X = np.random.multivariate_normal(mu, C, N)

# Normalize X (min-max scaling)
X = (X - np.min(X, axis=0)) / (np.max(X, axis=0) - np.min(X, axis=0))

# Define the function y = f(X)
y = 2 * X[:, 0] + X[:, 1] + 3 * np.exp(X[:, 0] * X[:, 1] * X[:, 2])

# Add random error
sigma2 = np.var(y) / 100  # Variance of random error
y += np.random.normal(0, np.sqrt(sigma2), N)

# HDMR options as a dictionary
options = {
    'graphics': 1,
    'maxorder': 3,
    'maxiter': 100,
    'bf1': 1,
    'bf2': 0,
    'bf3': 0,
    'M': 4,
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

	# Output for debugging or verification (displaying the first few rows of the result)
	print("Sensitivity results (S):\n", S)
