import numpy as np

# Convert correlation matrix to covariance matrix
def corr2cov(std_devs, corr_matrix):
    """
    Convert correlation matrix to covariance matrix.
    Assumes that 'std_devs' is a vector of standard deviations for each variable.
    """
    return np.outer(std_devs, std_devs) * corr_matrix
