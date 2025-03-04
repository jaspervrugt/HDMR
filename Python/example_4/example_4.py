% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE        4444   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE           44 44   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       44  44   %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE          444444   %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE          44   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Example from the following paper
%  Li, G., and H. Rabitz (2012), General formulation of HDMR component 
%      functions with independent and correlated variables, J. 
%      Math. Chem., 50, pp. 99-130

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
d = 3  # Number of parameters
N = 5000  # Number of samples

# Latin hypercube sampling
X = LH_sampling(np.zeros(3), np.ones(3), N)  # X samples from 0 to 1

# Normalize X values
X = (X - np.min(X, axis=0)) / (np.max(X, axis=0) - np.min(X, axis=0))

# Function to generate y = f(X)
y = np.sin(2 * np.pi * X[:, 0] - np.pi) + \
    7 * (np.sin(2 * np.pi * X[:, 1] - np.pi)) ** 2 + \
    0.1 * (2 * np.pi * X[:, 2] - np.pi) ** 4 * np.sin(2 * np.pi * X[:, 0] - np.pi)

# Variance of random error
sigma2 = np.var(y) / 55
y += np.random.normal(0, np.sqrt(sigma2), N)  # Add random error

# HDMR options (assuming you have a similar function to 'HDMR' in Python)
options = {
    'graphics': 1,
    'maxorder': 3,
    'maxiter': 100,
    'bf1': 1,
    'bf2': 0,
    'bf3': 0,
    'm': 4,
    'K': 100,
    'R': 500,
    'method': 1,
    'alfa': 0.01,
    'lambda': 0.10,
    'vartol': 1e-3,
    'refit': 1
}

if __name__ == '__main__':
	# Run the HDMR toolbox
	S, Ss, Fx, Em, Xy, RT = HDMR(X, y, options)

