import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from HDMR_functions import LH_sampling

def evaluate_TE(N, solver):
    # This function runs the Temperature - Extent model N different times

    # Specify temperature range, beta, initial state, and numerical solution
    plugin = {
        'T_e': np.arange(500, 701, 1),  # Temperature range from 500 to 700
        'beta': 1/6,                     # Constant beta
        'a0': 0,                         # Initial condition for state variable
        'solution': solver,              # Solver type: 'explicit' or 'implicit'
        'n': 201,                        # Number of temperature steps
        'options': {}                    # Options for ODE solver (empty here)
    }

    # Number of model parameters
    d = 2
    
    # Minimum and maximum values of parameters
    X_min = np.array([1.0, 1.00])
    X_max = np.array([500, 20.0])

    # Labels for parameters
    label_par = ['$E$', '$A$']

    # Sample randomly the parameter values within their ranges
    X = LH_sampling(X_min, X_max, N)

    # Initialize model output (N trials with n simulated counts)
    Y = np.full((N, plugin['n']), np.nan)

    # Initialize flags for checking valid simulations
    flag = np.full(N, np.nan)

    prev = 0
    # Now evaluate the model
    for i in range(N):
        _, a_s, flag[i] = Temp_extent(X[i, :2], plugin)
        # Store the model output
        Y[i, :] = a_s
        # Print progress
        prev = print(f'MODEL EVALUATIONS: {100 * (i / N):3.2f}% DONE', end='\r')

    # Remove "bad" simulations
    valid_idx = flag == 1
    Y = Y[valid_idx, :]
    X = X[valid_idx, :d]

    return X, label_par, Y

def Temp_extent(x, plugin):
    # Solves temperature - extent relationship for given E = x[0] and A = x[1]
    flag = 1  # "issues" flag
    R = 8.314e-3  # Universal gas constant

    # Two numerical solution methods for stiff differential equation
    if plugin['solution'].lower() == 'implicit':
        # Implicit solver using scipy's ode solver (ode23s equivalent)
        sol = solve_ivp(lambda t, a: ode_model(t, a, x, plugin['beta'], R), 
                        [plugin['T_e'][0], plugin['T_e'][-1]], [plugin['a0']], 
                        t_eval=plugin['T_e'], method='LSODA')
        T_s = sol.t
        a_s = sol.y[0, :]
    elif plugin['solution'].lower() == 'explicit':
        # Explicit solver: fixed dTee
        nT = 10000  # Number of temperature steps
        dT = (plugin['T_e'][-1] - plugin['T_e'][0]) / nT  # Integration step of temperature
        T = np.linspace(plugin['T_e'][0], plugin['T_e'][-1], nT + 1)
        T_s, a_s = explicit_model(T, plugin['T_e'], plugin['a0'], x, plugin['beta'], R, dT, nT)
    else:
        raise ValueError('Unknown integration method: {}'.format(plugin['solution']))

    # Now check output and set to zero if NaN
    if np.any(np.isnan(a_s)):
        a_s = np.zeros(plugin['n'])
        flag = -1

    return T_s, a_s, flag

def ode_model(t, a, x, beta, R):
    # Implicit solver function - returns da
    da = (1 / beta) * 10**x[1] * np.exp(-x[0] / R / t) * (1 - a)
    return da

def explicit_model(T, T_e, a0, x, beta, R, dT, nT):
    # Explicit solver integrates with dT
    a = np.full(nT + 1, np.nan)
    a[0] = a0  # Initial condition for state variable a

    for i in range(1, nT + 1):
        da = (1 / beta) * 10**x[1] * np.exp(-x[0] / R / T[i - 1]) * (1 - a[i - 1])
        a[i] = a[i - 1] + da * dT

    # Now interpolate [T, a] at T_e
    a_s = interp1d(T, a, kind='linear', fill_value="extrapolate")(T_e)
    return T, a_s
