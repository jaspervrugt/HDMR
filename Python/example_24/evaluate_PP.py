import numpy as np
from scipy.integrate import solve_ivp

def evaluate_PP(X, n, d, u0, ode_options, dt, t_max):
    """
    Evaluates the 1-predator-1-prey model.

    Parameters:
    X (np.array): Nxd matrix of parameter vectors
    n (int): Length of simulation (number of time steps)
    d (int): Number of parameters per vector (4 parameters in this case)
    u0 (np.array): Initial condition vector (prey and predator initial populations)
    ode_options (dict): Integration options (can be passed to solve_ivp)
    dt (float): Time step for the simulation
    t_max (float): Maximum time for the simulation

    Returns:
    Y (np.array): A 3D matrix of simulated prey and predator populations.
                  Shape: (N, n, K) where K = 2 (prey and predator).
    """

    N = X.shape[0]  # Number of parameter vectors
    K = 2            # Number of species (prey and predator)

    # Initialize output array: N simulations, each with n time steps, and K species
    Y = np.full((N, n, K), np.nan)
    
    prev = 0

    # Perform the simulation for each parameter vector in X
    for i in range(N):
        # Use solve_ivp to solve the system of ODEs (equivalent to ode45 in MATLAB)
        solution = solve_ivp(PP, [0, t_max], u0, t_eval=np.arange(0, t_max + dt, dt), args=(X[i, :d],), **ode_options)
        # Store the simulation results in the output array
        Y[i, :, :] = solution.y.T  # Transpose the result to match the shape (n, K)
        # Print progress
        prev = print(f'MODEL EVALUATIONS: {100 * (i + 1) / N:.2f}% DONE', end='\r')

    return Y


# Predator-Prey ODEs
def PP(u, t, params):
    # u[0] = prey population, u[1] = predator population
    r, alfa, m, theta = params
    K = 50  # carrying capacity for prey
    dxdt = r * u[0] * (1 - u[0] / K) - alfa * u[0] * u[1]
    dydt = -m * u[1] + theta * u[0] * u[1]
    return [dxdt, dydt]
