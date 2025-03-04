# ####################################################################### #
#                                                                         #
#   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE    11  11     #
#   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE        11  11     #
#   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE     11  11     #
#   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE        11  11     #
#   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE    11  11     #
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
from test_function import test_function

# Modified Sobol test function

#    Reference:
#    Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto,
#    M., and Tarantola, S. (2010). Variance-based sensitivity
#    analysis of model output. Design and estimator for the
#    total sensitivity index. Computer Physics Communications,
#    181, 259-270.

#    Configuration Used:
#    Campolongo, F., Saltelli, A., and Cariboni, J. (2011).
#    From screening to quantitative sensitivity analysis. A
#    unified approach. Computer Physics Communications,
#    182, 978-988.

d = 20  	# Number of parameters
N = 5000  	# Number of samples

# Sample X-values from U(0, 1)
X = np.random.rand(N, d)

# Initialize y values
y = np.empty(N)

# Loop through each sample to evaluate the test function
for i in range(N):
    y[i] = test_function(X[i])

# HDMR options dictionary (for use in your HDMR function)
options = {
    'graphics': 1,
    'maxorder': 3,
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

	# Print the results
	print(S, Ss, Fx, Em, Xy)

