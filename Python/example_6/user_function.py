import numpy as np

def user_function(X, mu, a, b, c, e):
    """
    Evaluate the user-defined function, y = f(X), for the Nxd matrix of parameter vectors X.
    
    Parameters:
    X : np.ndarray
        N x d matrix of input parameter vectors
    mu : np.ndarray
        Mean vector (mu) of length d
    a, b, c, e : np.ndarray
        Coefficients used in the function definitions
    
    Returns:
    y : np.ndarray
        Output vector calculated from the function
    """
    
    # Define the three functions g1, g2, and g3
    g1 = lambda x, mu, a, b: (a[1] * (x[:, 0] - mu[0]) + a[0]) * (b[1] * (x[:, 1] - mu[1]) + b[0])
    g2 = lambda x, mu, c: c[2] * (x[:, 1] - mu[1])**2 + c[1] * (x[:, 1] - mu[1]) + c[0]
    g3 = lambda x, mu, e: e[3] * (x[:, 2] - mu[2])**3 + e[2] * (x[:, 2] - mu[2])**2 + e[1] * (x[:, 2] - mu[2]) + e[0]
    
    # Compute the three components
    part1 = g1(X, mu, a, b)
    part2 = g2(X, mu, c)
    part3 = g3(X, mu, e)
    
    # The final result is the sum of the three parts
    y = part1 + part2 + part3
    
    return y


# Define the covariance matrix creation function
def create_covariance(s, r12):
    return np.array([[s[0]**2, r12*s[0]*s[1], 0],
                     [r12*s[0]*s[1], s[1]**2, 0],
                     [0, 0, s[2]**2]])
