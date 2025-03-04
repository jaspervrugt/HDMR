import numpy as np
from scipy.integrate import odeint
import multiprocessing
import os

def VDM_model(X, n, d, u0, ode_options, dt, t_end, mV, tbx):
    """
    Simulates the counts and summary metrics of the Vandermeer predator-prey model.

    Args:
        X (ndarray): Nxd matrix of parameter vectors.
        n (int): Length of simulation.
        d (int): Number of parameters.
        u0 (list): Initial conditions for 4 species.
        ode_options (dict): Options for ODE solver.
        dt (float): Time step size.
        t_end (float): End time for simulation.
        mV (int): Model version.
        tbx (int): Whether to run in parallel (1 for parallel, 0 for sequential).

    Returns:
        X (ndarray): The parameter vectors.
        Y (ndarray): The simulated counts.
        SS (ndarray): Summary statistics for each simulation.
    """
VDM(t, u_in, par, parameterization, tbx):
    N = X.shape[0]  # Number of parameter vectors
    Y = np.full((N, n, 4), np.nan)  # Initialize model output (N simulations, n time steps, 4 species)

    if tbx == 0:  # Sequential computation
        prev = 0
        for i in range(N):
            t = np.arange(0, t_end + dt, dt)  # Time array from 0 to t_end with step dt
            sol = odeint(VDM, u0, t, args=(X[i, :d], mV), **ode_options)  # Solve ODE system
            Y[i, :, :] = sol  # Store solution
            prev = print_progress(i, N, prev)
    
    elif tbx == 1:  # Parallel computation
        pool = multiprocessing.Pool()  # Initialize multiprocessing pool
        print(f"VDM: {pool._processes} workers are used to evaluate the model {N} times")

        def simulate(i):
            t = np.arange(0, t_end + dt, dt)
            sol = odeint(VDM, u0, t, args=(X[i, :d], mV), **ode_options)
            return sol

        results = pool.map(simulate, range(N))  # Parallelize the simulations
        for i in range(N):
            Y[i, :, :] = results[i]  # Store results
        pool.close()
        pool.join()

    # Summary statistics computation
    _, n, K = Y.shape
    SS = np.full((N, 8, K), np.nan)
    ct = 0

    # Loop over each species
    for k in range(K):
        for i in range(N):
            mxL, mxM = peakfinder(Y[i, :, k], 0.2)  # Find peaks
            if mxM.size > 0:
                # Calculate summary statistics for this simulation
                dim = mxM.shape
                mn_dis = n / (dim[1] + 1)  # Mean distance between peaks
                peak_mn = np.mean(mxM)  # Mean peak amplitude
                peak_mx = np.max(mxM)  # Max peak value
                
                SS[ct, :, k] = [
                    np.std(Y[i, :, k]),  # Standard deviation
                    np.mean(Y[i, :, k]),  # Mean value
                    len(mxL),  # Number of peaks
                    mn_dis,  # Mean distance between peaks
                    peak_mn,  # Mean peak amplitude
                    peak_mx,  # Max peak
                    X[i, 0],  # First parameter (beta)
                    X[i, 1]  # Second parameter (alpha)
                ]
                ct += 1

    return X, Y, SS


def VDM(t, u_in, par, parameterization):
    """
    Coupled ordinary differential equations for the Vandermeer predator-prey model.
    This function defines the system of differential equations for the four species.
    """
    # Assign the states
    P1, P2, Z1, Z2 = u_in

    if parameterization == 1:
        beta = par[0]
        alfa1 = par[1]
        alfa2 = alfa1
        r1 = par[2]
        r2 = par[3]
        g = par[4]
        K1 = par[5]
        K2 = K1
        m1 = par[6]
        m2 = m1
        H = par[7]
    elif parameterization == 2:
        beta = par[0]
        alfa1 = par[1]
        alfa2 = par[2]
        r1 = par[3]
        r2 = par[4]
        g = par[5]
        K1 = par[6]
        K2 = par[7]
        m1 = par[8]
        m2 = par[9]
        H = par[10]
    else:
        raise ValueError("Unknown parameterization")

    # Define the system of ODEs
    dP1dt = r1 * P1 * (1 - (1 / K1) * (P1 + alfa1 * P2)) - (
        g * P1 * Z1 / (H + (P1 + beta * P2)) + g * beta * P1 * Z2 / (H + (beta * P1 + P2))
    )
    dP2dt = r2 * P2 * (1 - (1 / K2) * (alfa2 * P1 + P2)) - (
        g * beta * P2 * Z1 / (H + (P1 + beta * P2)) + g * P2 * Z2 / (H + (beta * P1 + P2))
    )
    dZ1dt = g * (P1 + beta * P2) * Z1 / (H + (P1 + beta * P2)) - m1 * Z1
    dZ2dt = g * (beta * P1 + P2) * Z2 / (H + (beta * P1 + P2)) - m2 * Z2

    # Return the system of differential equations as a vector
    return [dP1dt, dP2dt, dZ1dt, dZ2dt]


def peakfinder(signal, threshold):
    """
    Function to find peaks in the signal that are above the given threshold.
    Returns the indices of the peaks (mxL) and the peak values (mxM).
    """
    mxL = []
    mxM = []

    for i in range(1, len(signal) - 1):
        if signal[i] > signal[i - 1] and signal[i] > signal[i + 1] and signal[i] > threshold:
            mxL.append(i)
            mxM.append(signal[i])

    return np.array(mxL), np.array(mxM)


def print_progress(i, N, prev):
    """
    Helper function to print the progress of the simulations.
    """
    prev = print(f'VDM: MODEL EVALUATIONS: {100*(i+1)/N:.2f}% DONE', end='\r')
    return prev
