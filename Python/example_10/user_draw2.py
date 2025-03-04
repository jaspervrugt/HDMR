import numpy as np
from scipy.stats import norm

# Define the user_draw2 function to generate samples from a normal mixture
def user_draw2(mu1, mu2, C1, C2, w1, N):
    """Generates samples from a normal mixture with two components."""
    # Draw samples for each normal distribution
    N1 = int(w1 * N)
    N2 = N - N1
    
    # Generate N1 samples from the first normal distribution (C1, mu1)
    X1 = np.random.multivariate_normal(mu1, C1, N1)
    
    # Generate N2 samples from the second normal distribution (C2, mu2)
    X2 = np.random.multivariate_normal(mu2, C2, N2)
    
    # Combine the two sets of samples
    X = np.vstack([X1, X2])
    
    return X