# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      999999   #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE          99  99   #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       999999   #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE              99   #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      999999   #
#                                                                         #
# ####################################################################### #

# Example 1 of the following paper
#  Chastaing, G., F. Gamboa, and C. Prieur (2012), Generalized Hoeffding-
#      Sobol decomposition for dependent variables - application to 
#      sensitivity analysis, Electron. J. Statistics, 6, pp. 2420-2448, 
#      https://doi.org/10.1214/12-EJS749

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
from sklearn.mixture import GaussianMixture
from user_draw import user_draw

# Case study
d = 2  					# number of dimensions
N = 1000  				# number of samples
mu1 = np.zeros((d, 1))  		# mean vector for the first normal distribution
mu2 = np.zeros((d, 1))  		# mean vector for the second normal distribution
w1 = 0.2  				# weight of the first normal distribution
r12 = 0.4  				# covariance between X1 and X2 for the second normal distribution
C1 = np.eye(d)  			# covariance matrix for the first normal distribution
C2 = np.array([[0.5, r12], [r12, 0.5]]) # covariance matrix for the second normal distribution

# Generate the samples X from the normal mixture
X = user_draw(mu1, mu2, C1, C2, w1, N, d)

# Compute the function y = f(X)
y = X[:, 0] + X[:, 1] + X[:, 0] * X[:, 1]

# Add random error to y
sigma2 = np.var(y) / 100  				# variance of the random error
y = y + np.random.normal(0, np.sqrt(sigma2), N)  	# add Gaussian noise

# Specify HDMR options (this part assumes you have a Python HDMR implementation)
options = {
    'graphics': 1,
    'maxorder': 2,
    'maxiter': 100,
    'bf1': 1,
    'bf2': 1,
    'bf3': 1,
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
