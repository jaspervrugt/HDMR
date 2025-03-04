import numpy as np

def B_spline(X, R, d, m):
    """
    Returns coefficients, B, of B-splines.
    
    Parameters:
    - X: RxD matrix, R vectors of d parameters
    - R: number of samples of X
    - d: number of parameters
    - m: number of B-spline intervals
    
    Returns:
    - B: Rx(m+3)xd matrix of B-spline coefficients
    """
    
    # Initialize matrix B
    B = np.zeros((R, m+3, d))
    
    # Calculate interval, h1
    h1 = 1 / m
    
    # Loop over each parameter of X and B-spline interval
    for i in range(d):
        for j in range(m+3):
            # Change value of k
            # k = j - 2
            k = j - 1
            for r in range(R):
                # Check interval of X[r,i]
                if X[r, i] > (k + 1) * h1 and X[r, i] <= (k + 2) * h1:
                    B[r, j, i] = ((k + 2) * h1 - X[r, i]) ** 3
                
                elif X[r, i] > k * h1 and X[r, i] <= (k + 1) * h1:
                    B[r, j, i] = ((k + 2) * h1 - X[r, i]) ** 3 - \
                                  4 * ((k + 1) * h1 - X[r, i]) ** 3
                
                elif X[r, i] > (k - 1) * h1 and X[r, i] <= k * h1:
                    B[r, j, i] = ((k + 2) * h1 - X[r, i]) ** 3 - \
                                  4 * ((k + 1) * h1 - X[r, i]) ** 3 + \
                                  6 * (k * h1 - X[r, i]) ** 3
                
                elif X[r, i] > (k - 2) * h1 and X[r, i] <= (k - 1) * h1:
                    B[r, j, i] = ((k + 2) * h1 - X[r, i]) ** 3 - \
                                  4 * ((k + 1) * h1 - X[r, i]) ** 3 + \
                                  6 * (k * h1 - X[r, i]) ** 3 - \
                                  4 * ((k - 1) * h1 - X[r, i]) ** 3
    
    # Multiply B with m^3
    B = m ** 3 * B
    
    return B
