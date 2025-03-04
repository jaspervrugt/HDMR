import numpy as np

def test_function(X):
    # Given parameters
    a = np.array([100, 0, 100, 100, 100, 100, 1, 10, 0, 0, 9, 0, 100, 100, 4, 100, 100, 7, 100, 2])
    alpha = np.array([1, 4, 1, 1, 1, 1, 0.4, 3, 0.8, 0.7, 2, 1.3, 1, 1, 0.3, 1, 1, 1.5, 1, 0.6])
    delt = np.array([0.2942, 0.2560, 0.3004, 0.5150, 0.7723, 0.4567, 0.8390, 0.1369,
                     0.1558, 0.4356, 0.0257, 0.3248, 0.0718, 0.9155, 0.6877, 0.5548, 
                     0.5835, 0.8083, 0.6309, 0.8071])

    # Initialize output list for each function value
    y = np.zeros(20)
    
    # Loop through the dimensions and calculate each y(i)
    for i in range(20):
        xi = X[i]  # X is expected to be a 1D array of length 20
        y[i] = ((1 + alpha[i]) * np.abs(2 * (xi + delt[i] - np.floor(xi + delt[i])) - 1) ** alpha[i] + a[i]) / (1 + a[i])
    
    # Return the product of all elements in y
    return np.prod(y)