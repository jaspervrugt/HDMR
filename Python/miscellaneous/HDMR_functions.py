# ####################################################################### #
#                                                                         #
#          HHH     HHH   DDDDDDDDD    MMM     MMM   RRRRRRRR              #
#          HHH     HHH   DDDDDDDDDD   MMM     MMM   RRRRRRRRR             #
#          HHH     HHH   DDD    DDD   MMMM   MMMM   RRR    RRR            #
#          HHH     HHH   DDD    DDD   MMMMM MMMMM   RRR    RRR            #
#          HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRRR             #
#          HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRR              #
#          HHH     HHH   DDD    DDD   MMM     MMM   RRRRRRRRR             #
#          HHH     HHH   DDD    DDD   MMM     MMM   RRR   RRR             #
#          HHH     HHH   DDDDDDDDDD   MMM     MMM   RRR    RRR            #
#          HHH     HHH   DDDDDDDDD    MMM     MMM   RRR     RRR           #
#                                                                         #
# ####################################################################### #
#                                                                         #
#  High-Dimensional Model Representation (HDMR) using B-spline functions  #
#  for variance-based global sensitivity analysis (GSA) with correlated   #
#  and uncorrelated inputs. This function uses as input a N x d matrix    #
#  of N different d-vectors of model inputs (factors/parameters) and a    #
#  N x 1 vector of corresponding model outputs and returns to the user    #
#  each factor's first, second, and third order sensitivity coefficient   #
#  (separated in total, structural and correlative contributions), an     #
#  estimate of their 95# confidence intervals (from bootstrap method)     #
#  and the coefficients of the significant B-spline basis functions that  #
#  govern output, y (determined by an F-test of the error residuals of    #
#  the HDMR model (emulator) with/without a given first, second and/or    #
#  third order B-spline). These coefficients define an emulator that can  #
#  be used to predict the output, y, of the original (CPU-intensive?)     #
#  model for any d-vector of model inputs. For uncorrelated model inputs  #
#  (columns of X are independent), the HDMR sensitivity indices reduce    #
#  to a single index (= structural contribution), consistent with their   #
#  values derived from commonly used variance-based GSA methods.          #
#                                                                         #
# ####################################################################### #
#                                                                         #
#  MAIN REFERENCE                                                         #
#   Li, G. H. Rabitz, P.E. Yelvington, O.O. Oluwole, F. Bacon,            #
#       C.E. Kolb, and J. Schoendorf (2010), Global sensitivity analysis  #
#       for systems with independent and/or correlated inputs, Journal of #
#       Physical Chemistry A, Vol. 114 (19), pp. 6022 - 6032, 2010        #
#   Gao, Y., A. Sahin, and J.A. Vrugt (2023), Probabilistic Sensitivity   #
#       Analysis With Dependent Variables: Covariance-Based Decomposition #
#       of Hydrologic Models, Water Resources Research, 59 (4),           #
#       e2022WR0328346, https://doi.org/10.1029/2022WR032834              #
#                                                                         #
#  PYTHON CODE                                                            #
#  © Written by Jasper A. Vrugt using GPT-4 OpenAI's language model       # 
#    University of California Irvine                                      #
#  Version 2.0    Dec 2024                                                #
#                                                                         #
# ####################################################################### #

import numpy as np                                      # type: ignore
import sys
import os
import itertools
from scipy.stats import f                               # type: ignore
from scipy.special import erfcinv                       # type: ignore
from B_spline import B_spline
import matplotlib.pyplot as plt                         # type: ignore
import pandas as pd                                     # type: ignore
from matplotlib.backends.backend_pdf import PdfPages    # type: ignore
from matplotlib.lines import Line2D                     # type: ignore
from screeninfo import get_monitors                     # type: ignore

def HDMR_setup(X, y, user_options):
    """
    Setup function for HDMR, including validation and configuration of options.

    Parameters:
    -----------
    X : numpy.ndarray
        Nxd matrix: N vectors of d parameters
    y : numpy.ndarray
        Nx1 vector: single model output for each row of X
    user_options : dict, optional
        User-defined options as a dictionary (if not provided, defaults will be used).

    Returns:
    --------
    N, d, graphics, maxorder, maxiter, bf1, bf2, bf3, m, K, R, method, alfa, vartol, lambda_, refit : various configuration parameters
    """

    # Check the shape of X
    N, d = X.shape
    if d == 1:
        raise ValueError("HDMR ERROR: Matrix X contains only a single column: No point to do sensitivity analysis when d = 1")
    if N < 100:
        raise ValueError(f"HDMR ERROR: Number of samples, N ({N}), is insufficient")
    
    # Check the shape of y
    if y.shape[0] != N:
        raise ValueError(f"HDMR ERROR: Dimension mismatch. The number of rows, N ({N}), of y should match number of rows of X")

    # Define default options
    def_options = {
        'graphics': 1,
        'maxorder': 3,
        'maxiter': 100,
        'bf1': 1,
        'bf2': 1,
        'bf3': 1,
        'm': 2,
        'K': 100,
        'R': N // 2,
        'method': 1,
        'alfa': 0.01,
        'lambda': 0.1,
        'vartol': 1e-5,
        'refit': 1
    }

    # Output configuration settings (equivalent to the file writing in MATLAB)
    with open('HDMR_settings.txt', 'w') as fid:
        fid.write('|-------------------------------------------------------------------------|\n')
        fid.write('|                                                                         |\n')
        fid.write('|           HHH     HHH   DDDDDDDDD    MMM     MMM   RRRRRRRR             |\n')
        fid.write('|           HHH     HHH   DDDDDDDDDD   MMM     MMM   RRRRRRRRR            |\n')
        fid.write('|           HHH     HHH   DDD    DDD   MMMM   MMMM   RRR    RRR           |\n')
        fid.write('|           HHH     HHH   DDD    DDD   MMMMM MMMMM   RRR    RRR           |\n')
        fid.write('|           HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRRR            |\n')
        fid.write('|           HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRR             |\n')
        fid.write('|           HHH     HHH   DDD    DDD   MMM     MMM   RRRRRRRRR            |\n')
        fid.write('|           HHH     HHH   DDD    DDD   MMM     MMM   RRR   RRR            |\n')
        fid.write('|           HHH     HHH   DDDDDDDDDD   MMM     MMM   RRR    RRR           |\n')
        fid.write('|           HHH     HHH   DDDDDDDDD    MMM     MMM   RRR     RRR          |\n')
        fid.write('|                                                                         |\n')
        fid.write('|                                                                         |\n')
        fid.write('|  SSSSSSS EEEEEEEE TTTTTTTT TTTTTTTT IIIIIIII NN    NN  GGGGGG   SSSSSSS |\n')
        fid.write('| SSSSSSS  EEEEEEEE TTTTTTTT TTTTTTTT  IIIIII  NNN   NN GGGGGG   SSSSSSS  |\n')
        fid.write('| SS       EE       TT TT TT TT TT TT    II    NNNN  NN GG       SS       |\n')
        fid.write('| SSSSSS   EEEEE    T  TT  T T  TT  T    II    NN NN NN GG  GGGG SSSSSS   |\n')
        fid.write('| SSSSSSS  EEEEE       TT       TT       II    NN NN NN GG   GGG SSSSSSS  |\n')
        fid.write('|       SS EE          TT       TT       II    NN  NNNN GG    GG       SS |\n')
        fid.write('|  SSSSSSS EEEEEEEE    TT       TT     IIIIII  NN   NNN GGGGGGGG  SSSSSSS |\n')
        fid.write('| SSSSSSS  EEEEEEEE    TT       TT    IIIIIIII NN    NN  GGGGGG  SSSSSSS  |\n')
        fid.write('|                                                                         |\n')
        fid.write('|-------------------------------------------------------------------------|\n')
    
        # Ensure that user options are within valid bounds
        graphics = int(user_options.get('graphics', def_options['graphics']))
        if graphics not in [0, 1]:
            raise ValueError('HDMR_EXT ERROR: Field "graphics" of options should take on the value of 0 or 1')

        maxorder = int(user_options.get('maxorder', def_options['maxorder']))
        if maxorder not in [1, 2, 3]:
            raise ValueError('HDMR ERROR: Field "maxorder" of options should be an integer with values of 1, 2 or 3')
    
        if d == 2 and maxorder > 2:
            print("HDMR WARNING: Field 'maxorder' of options set to 2 as d = 2 (X has two columns)")
            maxorder = 2

        maxiter = int(user_options.get('maxiter', def_options['maxiter']))
        if not (1 <= maxiter <= 1000):
            raise ValueError('HDMR ERROR: Field "maxiter" of options should be an integer between 1 and 1000')

        bf1 = int(user_options.get('bf1', def_options['bf1']))
        if bf1 not in [0, 1]:
            raise ValueError('HDMR ERROR: Field "bf1" of options should be an integer with value 0 or 1')

        bf2 = int(user_options.get('bf2', def_options['bf2']))
        if bf2 not in [0, 1]:
            raise ValueError('HDMR ERROR: Field "bf2" of options should be an integer with value 0 or 1')

        bf3 = int(user_options.get('bf3', def_options['bf3']))
        if bf3 not in [0, 1]:
            raise ValueError('HDMR ERROR: Field "bf3" of options should be an integer with value 0 or 1')

        m = int(user_options.get('m', def_options['m']))
        if not (1 <= m <= 10):
            raise ValueError('HDMR ERROR: Field "m" of options should be an integer between 1 and 10')

        K = int(user_options.get('K', def_options['K']))
        if not (1 <= K <= 500):
            raise ValueError('HDMR ERROR: Field "K" of options should be an integer between 1 and 500')

        R = int(user_options.get('R', def_options['R']))
        if not (100 <= R <= N):
            raise ValueError(f'HDMR ERROR: Field "R" of options should be an integer between 100 and N ({N}), the number of rows matrix X')

        method = int(user_options.get('method', def_options['method']))
        if method not in [1, 2]:
            raise ValueError('HDMR ERROR: Field "method" of options should take on the value of 1 (forward selection) or 2 (backward elimination)')

        alfa = user_options.get('alfa', def_options['alfa'])
        if isinstance(alfa, str):
            raise ValueError('HDMR ERROR: Field "alfa" (confidence interval) of options should not be a string but a numerical value')
        if alfa > 0.1 or alfa < np.finfo(float).eps:
            raise ValueError('HDMR ERROR: Field "alfa" (confidence interval) of options should not exceed 0.1 or be smaller than eps')

        lambda_ = user_options.get('lambda', def_options['lambda'])
        if isinstance(lambda_, str):
            raise ValueError('HDMR ERROR: Field "lambda" (regularization term) of options should not be a string but a numerical value')
        if lambda_ > 10:
            print("HDMR WARNING: Field 'lambda' of options set rather large. Default: lambda = 0.1")
            lambda_ = 0.1
        if lambda_ < 0:
            raise ValueError('HDMR ERROR: Field "lambda" (regularization term) of options cannot be smaller than zero')

        vartol = user_options.get('vartol', def_options['vartol'])
        if isinstance(vartol, str):
            raise ValueError('HDMR ERROR: Field "vartol" (convergence threshold) of options should not be a string but a numerical value')
        if not (0 < vartol < 0.1):
            raise ValueError('HDMR ERROR: Field "vartol" (convergence threshold) of options should be between 0 and 0.1')

        refit = int(user_options.get('refit', def_options['refit']))
        if refit not in [0, 1]:
            raise ValueError('HDMR ERROR: Field "refit" of options should take on the value of 0 or 1')

        # Save options to file
        fid.write(f'\n          |===================================|\n')
        fid.write(f'          |  field of options      value      |\n')
        fid.write(f'          |-----------------------------------|\n')
        fid.write(f'          |    graphics   \t {graphics:8d}     |\n')
        fid.write(f'          |    maxorder   \t {maxorder:8d}     |\n')
        fid.write(f'          |    maxiter    \t {maxiter:8d}     |\n')
        fid.write(f'          |    bf1        \t {bf1:8d}     |\n')
        fid.write(f'          |    bf2        \t {bf2:8d}     |\n')
        fid.write(f'          |    bf3        \t {bf3:8d}     |\n')
        fid.write(f'          |    method     \t {method:8d}     |\n')       
        fid.write(f'          |    m          \t {m:8d}     |\n')
        fid.write(f'          |    K          \t {K:8d}     |\n')
        fid.write(f'          |    R          \t {R:8d}     |\n')
        fid.write(f'          |    alfa       \t {alfa:8.2f}     |\n')
        fid.write(f'          |    lambda     \t {lambda_:8f}     |\n')
        fid.write(f'          |    vartol     \t {vartol:8.2e}     |\n')
        fid.write(f'          |    refit      \t {refit:8d}     |\n')
        fid.write(f'          |===================================|\n')

    # Return the parameters
    return N, d, graphics, maxorder, maxiter, bf1, bf2, bf3, m, K, R, method, alfa, vartol, lambda_, refit


def HDMR_initialize(X, y, N, d, K, R, m, maxorder):
    # Set random seed
    np.random.seed(1 + round(100 * np.random.rand()))

    # STRUCTURE XY: Define content
    id = np.argsort(np.random.rand(N, K), axis=0)
    Xy = {
        'X_n': np.nan * np.ones((N, d)),
        'minX': np.min(X, axis=0),
        'maxX': np.max(X, axis=0),
        'y': y,
        'R': R,
        'id': id[:R, :K]
        }

    # Compute normalized X-values
    Xy['X_n'] = (X - Xy['minX']) / (Xy['maxX'] - Xy['minX'])

    # STRUCTURE Em: Important variables
    n2, n3 = 0, 0
    c2, c3 = [], []
    c1 = np.arange(1, d+1)
    n1 = d

    if maxorder > 1:
        c2 = np.array(list(itertools.combinations(range(0, d), 2)))
        n2 = len(c2)
    if maxorder == 3:
        c3 = np.array(list(itertools.combinations(range(0, d), 3)))
        n3 = len(c3)
    
    # Calculate total number of coefficients
    n = n1 + n2 + n3

    # Initialize m1, m2, and m3
    m1 = m + 3
    m2 = m1**2
    m3 = m1**3

    # STRUCTURE Em: Initialization
    Em = {}

    if maxorder == 1:
        Em = {
            'nterms': np.nan * np.ones(K),
            'p0': np.nan * np.ones(K),
            'RMSE': np.nan * np.ones(K),
            'm': m,
            'Y_e': np.nan * np.ones((R, K)),
            'f0': np.nan * np.ones(K),
            'c1': c1,
            'n1': d,
            'c2': c2,
            'n2': n2,
            'c3': c3,
            'n3': n3,
            'n': n,
            'maxorder': maxorder,
            'select': np.nan * np.ones((n, K)),
            'C1': np.nan * np.ones((m1, n1, K)),
            'C2': np.nan * np.ones((1, 1, K)),
            'C3': np.nan * np.ones((1, 1, K)),
            'B1': np.zeros((N, m1, n1)),
            'B2': np.nan * np.ones((N, 1)),
            'B3': np.nan * np.ones((N, 1)),
            'iter': np.nan * np.ones((4, K))
        }
    elif maxorder == 2:
        Em = {
            'nterms': np.nan * np.ones(K),
            'p0': np.nan * np.ones(K),
            'RMSE': np.nan * np.ones(K),
            'm': m,
            'Y_e': np.nan * np.ones((R, K)),
            'f0': np.nan * np.ones(K),
            'c1': c1,
            'n1': n1,
            'c2': c2,
            'n2': n2,
            'c3': c3,
            'n3': n3,
            'n': n,
            'maxorder': maxorder,
            'select': np.nan * np.ones((n, K)),
            'C1': np.nan * np.ones((m1, n1, K)),
            'C2': np.nan * np.ones((m2, n2, K)),
            'C3': np.nan * np.ones((1, 1, K)),
            'B1': np.zeros((N, m1, n1)),
            'B2': np.zeros((N, m2, n2)),
            'B3': np.nan * np.ones((N, 1)),
            'iter': np.nan * np.ones((4, K))
        }
    elif maxorder == 3:
        Em = {
            'nterms': np.nan * np.ones(K),
            'p0': np.nan * np.ones(K),
            'RMSE': np.nan * np.ones(K),
            'm': m,
            'Y_e': np.nan * np.ones((R, K)),
            'f0': np.nan * np.ones(K),
            'c1': c1,
            'n1': n1,
            'c2': c2,
            'n2': n2,
            'c3': c3,
            'n3': n3,
            'n': n,
            'maxorder': maxorder,
            'select': np.nan * np.ones((n, K)),
            'C1': np.nan * np.ones((m1, n1, K)),
            'C2': np.nan * np.ones((m2, n2, K)),
            'C3': np.nan * np.ones((m3, n3, K)),
            'B1': np.zeros((N, m1, n1)),
            'B2': np.zeros((N, m2, n2)),
            'B3': np.zeros((N, m3, n3)),
            'iter': np.nan * np.ones((4, K))
        }

    # Now compute B-spline values for all N samples of X_n
    Em['B1'] = B_spline(Xy['X_n'], N, d, m)

    # Now compute B values for second order
    if maxorder > 1:
        beta = [p for p in itertools.product(range(0, m1), repeat=2)]
        for k in range(n2):
            for j in range(m2):
                Em['B2'][0:N, j, k] = Em['B1'][0:N, beta[j][0], Em['c2'][k][0]] * \
                                      Em['B1'][0:N, beta[j][1], Em['c2'][k][1]]

    # Compute B values for third order
    if maxorder == 3:
        beta = [p for p in itertools.product(range(0, m1), repeat=3)]
        for k in range(n3):
            for j in range(m3):
                Em['B3'][0:N, j, k] = Em['B1'][0:N, beta[j][0], Em['c3'][k][0]] * \
                                    Em['B1'][0:N, beta[j][1], Em['c3'][k][1]] * \
                                    Em['B1'][0:N, beta[j][2], Em['c3'][k][2]]

    # STRUCTURE SA: Sensitivity analysis and analysis of variance decomposition
    SA = {
        'S': np.nan * np.ones((Em['n'], K)),
        'Sa': np.nan * np.ones((Em['n'], K)),
        'Sb': np.nan * np.ones((Em['n'], K)),
        'ST': np.nan * np.ones((d, 1)),
        'V_em': np.nan * np.ones((Em['n'], K)),
        'V_y': np.nan * np.ones(K)
    }

    # Return runtime
    RT = np.zeros(K)
    # Initialize various terms
    Y_em = np.nan * np.ones((R, Em['n']))
    T2, T3 = [], []
    
    # Initialize index of first, second, and third order terms (columns of Y_bf)
    j1 = np.arange(0, n1)
    j2 = np.arange(n1, n1 + n2)
    j3 = np.arange(n1 + n2, n)
    
    # Now initialize number of iterations first, second, and third order
    it2, it3, itr = 0, 0, 0

    # Print to screen
    print('  ----------------------------------------------------             ')
    print("  HHH     HHH   DDDDDDDDD    MMM     MMM   RRRRRRRR                ")
    print("  HHH     HHH   DDDDDDDDDD   MMM     MMM   RRRRRRRRR               ")
    print("  HHH     HHH   DDD    DDD   MMMM   MMMM   RRR    RRR              ")
    print("  HHH     HHH   DDD    DDD   MMMMM MMMMM   RRR    RRR              ")
    print("  HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRRR       /^ ^\   ")
    print("  HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRR       / 0 0 \  ")
    print("  HHH     HHH   DDD    DDD   MMM     MMM   RRRRRRRRR      V\ Y /V  ")
    print("  HHH     HHH   DDD    DDD   MMM     MMM   RRR   RRR       / - \   ")
    print("  HHH     HHH   DDDDDDDDDD   MMM     MMM   RRR    RRR     /     |  ")
    print("  HHH     HHH   DDDDDDDDD    MMM     MMM   RRR     RRR    V__) ||  ")
    print('  ----------------------------------------------------             ')
    print('  © Jasper A. Vrugt, University of California Irvine & GPT-4 OpenAI''s language model')
    print('    ________________________________________________________________________')
    print('    Version 2.0, Dec. 2024, Beta-release: MATLAB implementation is benchmark')
    print('\n')

    return Xy, Em, SA, RT, Y_em, T2, T3, m1, m2, m3, j1, j2, j3, it2, it3, itr


def HDMR_1st(B1, y_res, R, n1, m1, vartol, lambda_val, maxiter, bf):
    # Initialize coefficients, first-order contributions, temporary matrix, and iteration counter
    C1 = np.zeros((m1, n1))     # Coefficients
    y_i = np.zeros((R, n1))     # First-order contributions
    T1 = np.zeros((m1, R, n1))  # Temporary matrix for first-order terms
    it = 0  # Iteration counter

    # First-order individual estimation
    for j in range(n1):
        # Regularized least squares inversion (minimize || C1 ||_2)
        B11 = np.dot(B1[:R, :m1, j].T, B1[:R, :m1, j])
        cond_number = np.linalg.cond(B11)
        if cond_number > 1e15:
            pass
        else:
            T1[:m1, :R, j] = np.linalg.solve(B11 + lambda_val * np.eye(m1), B1[:R, :m1, j].T)
        C1[:m1, j] = np.dot(T1[:m1, :R, j], y_res)
        y_i[:R, j] = np.dot(B1[:R, :m1, j], C1[:m1, j])

    # First-order backfitting if bf == 1
    if bf == 1:
        var1b_old = np.sum(C1**2, axis=0)  # Initial variance
        varmax = 1

        while varmax > vartol and it < maxiter:
            for j in range(n1):
                y_r = y_res.copy()

                # Remove contributions from all other first-order terms
                for z in range(n1):
                    if j != z:
                        y_r -= np.dot(B1[:R, :m1, z], C1[:m1, z])

                # Regularized least squares inversion (minimize || C1 ||_2)
                C1[:m1, j] = np.dot(T1[:m1, :R, j], y_r)

            var1b_new = np.sum(C1**2, axis=0)
            varmax = np.max(np.abs(var1b_new - var1b_old))  # Check for convergence
            var1b_old = var1b_new
            it += 1

        # Now compute first-order terms
        for j in range(n1):
            y_i[:R, j] = np.dot(B1[:R, :m1, j], C1[:m1, j])
            # Subtract each first-order term from residuals
            y_res -= y_i[:R, j]

    return y_res, y_i, C1, T1, it


def HDMR_2nd(B2, y_res, R, n2, m2, vartol, lambda_val, maxiter, bf):
    # Initialize coefficients, second-order contributions, temporary matrix, and iteration counter
    C2 = np.zeros((m2, n2))     # Coefficients
    y_ij = np.zeros((R, n2))    # Second-order contributions
    T2 = np.zeros((m2, R, n2))  # Temporary matrix for second-order terms
    it = 0  # Iteration counter

    # Second-order individual estimation
    for j in range(n2):
        # Regularized least squares inversion (minimize || C2 ||_2)
        B22 = np.dot(B2[:R, :m2, j].T, B2[:R, :m2, j])
        cond_number = np.linalg.cond(B22)
        if cond_number > 1e15:
            pass
        else:
            T2[:m2, :R, j] = np.linalg.solve(B22 + lambda_val * np.eye(m2), B2[:R, :m2, j].T)
        C2[:m2, j] = np.dot(T2[:m2, :R, j], y_res)
        y_ij[:R, j] = np.dot(B2[:R, :m2, j], C2[:m2, j])

    # Second-order backfitting if bf == 1
    if bf == 1:
        var2b_old = np.sum(C2**2, axis=0)  # Initial variance
        varmax = 1

        while varmax > vartol and it < maxiter:
            for j in range(n2):
                y_r = y_res.copy()

                # Remove contributions from all other second-order terms
                for z in range(n2):
                    if j != z:
                        y_r -= np.dot(B2[:R, :m2, z], C2[:m2, z])

                # Regularized least squares inversion (minimize || C2 ||_2)
                C2[:m2, j] = np.dot(T2[:m2, :R, j], y_r)

            var2b_new = np.sum(C2**2, axis=0)
            varmax = np.max(np.abs(var2b_new - var2b_old))  # Check for convergence
            var2b_old = var2b_new
            it += 1

        # Now compute second-order terms
        for j in range(n2):
            y_ij[:R, j] = np.dot(B2[:R, :m2, j], C2[:m2, j])
            # Subtract each second-order term from residuals
            y_res -= y_ij[:R, j]

    return y_res, y_ij, C2, T2, it

def HDMR_3rd(B3, y_res, R, n3, m3, vartol, lambda_val, maxiter, bf):
    # Initialize coefficients, third-order contributions, temporary matrix, and iteration counter
    C3 = np.zeros((m3, n3))     # Coefficients
    y_ijk = np.zeros((R, n3))   # Third-order contributions
    T3 = np.zeros((m3, R, n3))  # Temporary matrix for third-order terms
    it = 0                      # Iteration counter

    # Third-order individual estimation
    for j in range(n3):
        # Regularized least squares inversion (minimize || C3 ||_2)
        B33 = np.dot(B3[:R, :m3, j].T, B3[:R, :m3, j])
        cond_number = np.linalg.cond(B33)
        if cond_number > 1e15:
            pass
        else:
            T3[:m3, :R, j] = np.linalg.solve(B33 + lambda_val * np.eye(m3), B3[:R, :m3, j].T)
        C3[:m3, j] = np.dot(T3[:m3, :R, j], y_res)
        y_ijk[:R, j] = np.dot(B3[:R, :m3, j], C3[:m3, j])

    # Third-order backfitting if bf == 1
    if bf == 1:
        var3b_old = np.sum(C3**2, axis=0)  # Initial variance
        varmax = 1

        while varmax > vartol and it < maxiter:
            for j in range(n3):
                y_r = y_res.copy()

                # Remove contributions from all other third-order terms
                for z in range(n3):
                    if j != z:
                        y_r -= np.dot(B3[:R, :m3, z], C3[:m3, z])

                # Regularized least squares inversion (minimize || C3 ||_2)
                C3[:m3, j] = np.dot(T3[:m3, :R, j], y_r)

            var3b_new = np.sum(C3**2, axis=0)
            varmax = np.max(np.abs(var3b_new - var3b_old))  # Check for convergence
            var3b_old = var3b_new
            it += 1

        # Now compute third-order terms
        for j in range(n3):
            y_ijk[:R, j] = np.dot(B3[:R, :m3, j], C3[:m3, j])
            # Subtract each third-order term from residuals (not needed in this case)
            # y_res -= y_ijk[:R, j]

    return y_ijk, C3, T3, it


def HDMR_F_test(Y, f0, Y_bf, R, alfa, m1, m2, m3, n1, n2, n3, n, method):
    # Initialize ind with zeros (all terms insignificant)
    ind = np.zeros(n)
    
    # Determine the significant components of the HDMR model via the F-test
    if method == 1:  # forward selection of terms (start with Y = f0)
        Y_res0 = Y - f0
        SSR0 = np.sum(Y_res0**2)
        p0 = 0

        for i in range(n):
            # model with ith term included
            Y_res1 = Y_res0 - Y_bf[:, i]
            
            # Number of parameters of proposed model (order dependent)
            if i < n1:
                p1 = p0 + m1  # 1st order
            elif i >= n1 and i < n1 + n2:
                p1 = p0 + m2  # 2nd order
            else:
                p1 = p0 + m3  # 3rd order

            # Calculate SSR of Y1
            SSR1 = np.sum(Y_res1**2)

            # Now calculate the F_stat (F_stat > 0 -> SSR1 < SSR0 )
            F_stat = ((SSR0 - SSR1) / (p1 - p0)) / (SSR1 / (R - p1))

            # Now calculate critical F value at confidence level 1 - alfa
            F_crit = f.ppf(1 - alfa, p1 - p0, R - p1)

            # Now determine whether to accept ith component into model
            if F_stat > F_crit:
                # ith term is significant and should be included in model
                ind[i] = 1
                Y_res0 = Y_res1
                SSR0 = SSR1
                p0 = p1

    elif method == 2:  # backward elimination of terms (start with Y = f0 + sum(all_terms))
        # NOTE: ONLY POSSIBLE IF R - p1 > 0 OTHERWISE MUST DO FORWARD SELECTION
        # Determine output of full model (all terms included)
        Y_res0 = Y - f0 - np.sum(Y_bf, axis=1)
        
        # Calculate the SSR of the full model (all terms included)
        SSR0 = np.sum(Y_res0**2)
        
        # Determine number of parameters of full model
        p0 = n1 * m1 + n2 * m2 + n3 * m3
        
        for i in range(n):
            # previous model with ith term excluded
            Y_res1 = Y_res0 + Y_bf[:, i]

            # Number of parameters of proposed model (order dependent)
            if i < n1:
                p1 = p0 - m1  # 1st order
            elif i >= n1 and i < n1 + n2:
                p1 = p0 - m2  # 2nd order
            else:
                p1 = p0 - m3  # 3rd order

            # Calculate SSR of Y1
            SSR1 = np.sum(Y_res1**2)

            # Now calculate the F_stat (F_stat > 0 if SSR1 > SSR0 as p1 < p0)
            F_stat = ((SSR0 - SSR1) / (p1 - p0)) / (SSR1 / (R - p1))

            # Now calculate critical F value at confidence level 1 - alfa
            F_crit = f.ppf(1 - alfa, p0 - p1, R - p1)  # Note 2nd term turned around

            # Now determine whether ith component is significant
            if F_stat > F_crit:
                # ith term is significant and should retain in model
                ind[i] = 1
            else:
                # ith term is insignificant and will be removed from model
                p0 = p1
                SSR0 = SSR1
                Y_res0 = Y_res1

    # Now return number of terms of final HDMR model
    nterms = np.sum(ind)
    
    return ind, nterms, p0


def HDMR_refit(B1, B2, B3, T1, T2, T3, C1, C2, C3, y_res, Y_em, R, n1, n2, n, m1, m2, m3, j1, j2, j3, ind, vartol, maxiter, maxorder):
    # Initialize important variables used for convergence analysis backfitting
    varb_old = np.nan * np.ones(n)
    varb_max = 1
    iter = 0
    i2 = i3 = []

    # Set insignificant 1st order terms to zero
    i1 = np.where(ind[j1] == 0)[0]
    C1[0:m1, i1] = 0
    Y_em[0:R, j1[i1]] = 0
    varb_old[j1] = np.sum(C1**2, axis=0)

    if maxorder > 1:
        # Set insignificant 2nd order terms to zero
        i2 = np.where(ind[j2] == 0)[0]
        C2[0:m2, i2] = 0
        Y_em[0:R, j2[i2]] = 0
        varb_old[j2] = np.sum(C2**2, axis=0)

    if maxorder > 2:
        # Set insignificant 3rd order terms to zero
        i3 = np.where(ind[j3] == 0)[0]
        C3[0:m3, i3] = 0
        Y_em[0:R, j3[i3]] = 0
        varb_old[j3] = np.sum(C3**2, axis=0)

    # Collect all coefficients with insignificant sensitivity
    ii = np.concatenate([i1, j2[i2], j3[i3]])

    # Determine index of all significant terms
    id = np.setdiff1d(np.arange(n), ii)

    # Backfitting of all orders combined
    while varb_max > vartol and iter < maxiter:
        varb_new = np.zeros(n)
        
        for i in range(len(id)):
            # Get column and index of all other non-zero columns
            j = id[i]
            ii = id.copy()
            ii = np.delete(ii, i)
            
            # Now remove all non-zero columns (nz) from Y_em except column j
            y_nzl = y_res - np.sum(Y_em[:, ii], axis=1)
            
            if j in j1:
                C1[0:m1, j] = T1[0:m1, 0:R, j] @ y_nzl[0:R]
                Y_em[0:R, j] = B1[0:R, 0:m1, j] @ C1[0:m1, j]
                varb_new[j] = np.sum(C1[0:m1, j]**2)
            elif j in j2:
                z = j - n1
                C2[0:m2, z] = T2[0:m2, 0:R, z] @ y_nzl[0:R]
                Y_em[0:R, j] = B2[0:R, 0:m2, z] @ C2[0:m2, z]
                varb_new[j] = np.sum(C2[0:m2, z]**2)
            else:
                z = j - n1 - n2
                C3[0:m3, z] = T3[0:m3, 0:R, z] @ y_nzl[0:R]
                Y_em[0:R, j] = B3[0:R, 0:m3, z] @ C3[0:m3, z]
                varb_new[j] = np.sum(C3[0:m3, z]**2)

        varb_max = np.max(np.abs(varb_new - varb_old))
        varb_old = varb_new
        iter += 1

    return Y_em, C1, C2, C3, iter


def ANCOVA(y, Y_em, V_y, R, n):
    """
    Analysis of Covariance (ANCOVA)
    Computes sensitivity indices for each emulator term.
    
    Parameters:
    - y: Rx1 vector, single model output for each row of X
    - Y_em: Rxn matrix with HDMR model terms from backfitting
    - V_y: scalar, variance of the model output
    - R: number of samples of X
    - n: number of terms of HDMR model
    
    Returns:
    - S: nx1 vector of total sensitivity for each emulator term
    - S_a: nx1 vector of structural sensitivity for each emulator term
    - S_b: nx1 vector of correlative sensitivity for each emulator term
    - V_em: nx1 vector of variance for each emulator term
    """
    
    # Compute the sum of all Y_bf terms
    Y0 = np.sum(Y_em, axis=1)  # sum of all Y_em terms across columns for each row
    # Initialize each variable
    S = np.full((n, 1), np.nan)
    S_a = np.full((n, 1), np.nan)
    S_b = np.full((n, 1), np.nan)
    V_em = np.full((n, 1), np.nan)
    
    # Analysis of covariance -> extension of analysis of variance
    for j in range(n):
        # Covariance matrix of j-th term of Y_em and actual Y
        C = np.cov(Y_em[:, j], y, rowvar=False)
        S[j, 0] = C[0, 1] / V_y     # Total sensitivity of j-th term
        # Covariance matrix of j-th term with emulator Y without j-th term
        C = np.cov(Y_em[:, j], Y0 - Y_em[:, j], rowvar=False)
        S_a[j, 0] = C[0, 0] / V_y   # Structural contribution of j-th term
        S_b[j, 0] = C[0, 1] / V_y   # Correlative contribution of j-th term
        # Variance in Y of j-th term
        V_em[j, 0] = np.sum(C[0, 0:2])  # Variance of the j-th term

    S = S.reshape(-1)
    S_a = S_a.reshape(-1)
    S_b = S_b.reshape(-1)
    V_em = V_em.reshape(-1)

    return S, S_a, S_b, V_em


def HDMR_end(SA, nterms, p0, RMSE, select, K, C2, C3, n1, n2, n3, n):
    """
    Prepares the return arguments of HDMR (High Dimensional Model Representation)
    
    Parameters:
    - SA: Sensitivity analysis data matrix
    - nterms: Number of terms
    - p0: Initial point
    - RMSE: Root mean square error
    - select: Selection matrix
    - K: Number of bootstrap trials
    - C2, C3: Confidence intervals (not directly used in this function)
    - n1, n2, n3: Dimensions
    - n: The number of variables
    
    Returns:
    - SA: Sensitivity analysis matrix
    - SA_sig: Sensitivity analysis with significant terms
    - Fx: A matrix containing some results related to the emulator
    """
    # Set the alpha (significance level) for bootstrap confidence intervals
    alfa = 0.05

    # Step 1: Compile table, SI, with all sensitivity results
    f_ret = np.sum(select, axis=1)
    SA, Fx = construct_tables(SA, nterms, p0, RMSE, alfa, f_ret, K, C2, C3, n1, n2, n3, n)

    # Identify significant terms
    nr = len(SA)
    id_sig = np.concatenate(([1], np.where(f_ret > 0)[0] + 2, [nr])) - 1
    SA_sig = [SA[i][:13] for i in id_sig]

    # Step 2: Print a table with emulator results
    with open('HDMR_results.txt', 'w+') as fid:
        fid.write('|-------------------------------------------------------------------------|\n')
        fid.write('|                                                                         |\n')
        fid.write('|           HHH     HHH   DDDDDDDDD    MMM     MMM   RRRRRRRR             |\n')
        fid.write('|           HHH     HHH   DDDDDDDDDD   MMM     MMM   RRRRRRRRR            |\n')
        fid.write('|           HHH     HHH   DDD    DDD   MMMM   MMMM   RRR    RRR           |\n')
        fid.write('|           HHH     HHH   DDD    DDD   MMMMM MMMMM   RRR    RRR           |\n')
        fid.write('|           HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRRR            |\n')
        fid.write('|           HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRR             |          /^ ^\  \n')
        fid.write('|           HHH     HHH   DDD    DDD   MMM     MMM   RRRRRRRRR            |         / 0 0 \ \n')
        fid.write('|           HHH     HHH   DDD    DDD   MMM     MMM   RRR   RRR            |         V\ Y /V \n')
        fid.write('|           HHH     HHH   DDDDDDDDDD   MMM     MMM   RRR    RRR           |          / - \  \n')
        fid.write('|           HHH     HHH   DDDDDDDDD    MMM     MMM   RRR     RRR          |         /     | \n')
        fid.write('|                                                                         |        V__) ||  \n')
        fid.write('|-------------------------------------------------------------------------|\n')
        fid.write('\n')

        # Table header based on number of bootstrap trials
        if K == 1:
            fid.write(f'Table 1: Properties of HDMR emulator, y = f(x), for training data set (no bootstrapping)\n')
        else:
            fid.write(f'Table 1: Properties of HDMR emulator, y = f(x), for randomized training data set of {K} bootstrap trials\n')
        fid.write('=========================================\n')
        fid.write('Emulator    # terms   # coefs.   RMSE    \n')
        fid.write('-----------------------------------------\n')

        # Print emulator results
        fmt_1 = '   {:<7d}\t{:<3d}      {:<3d}     {:<5.3f}\n'
        for k in range(1, K + 1):
            # Assuming Fx is a list of lists, print based on the structure
            fid.write(fmt_1.format(int(Fx[k][0]), Fx[k][1], Fx[k][2], float(Fx[k][3])))
        fid.write('=========================================\n')
        fid.write('\n')

        # Print sensitivity results
        for pr_tab in range(1, 3):
            fid.write('\n\n')
            if pr_tab == 1:
                id_table = np.arange(nr)
                SAtable = SA
            else:
                id_table = id_sig
                nr = len(id_sig)
                SAtable = SA_sig

            # Table headers for sensitivity results
            if K == 1:
                if pr_tab == 1:
                    fid.write(f'Table 2: HDMR results for all model components using all X and Y data (no bootstrapping)\n')
                else:
                    fid.write(f'Table 3: HDMR results for significant model components only using all X and Y data (no bootstrapping)\n')
            elif K > 1:
                if pr_tab == 1:
                    fid.write(f'Table 2: HDMR results for all model components using {K} bootstrap trials\n')
                else:
                    fid.write(f'Table 3: HDMR results for significant model components only using {K} bootstrap trials\n')

            fid.write('============================================================================================================ \n')
            fid.write('                                         ( = JOINT DETERMINATION )                                    \n')
            fid.write('                  ------------------------------------------------------------------------------------------ \n')
            fid.write('Term       Order        S^a             S^b              S               ST          V(i)/V(Y)     #select   \n')
            fid.write('------------------------------------------------------------------------------------------------------------ \n')

            if K == 1:
                fmt_1 = '{:<11}  {:>1}     {:>6.3f} (-----)  {:>6.3f} (-----)  {:>6.3f} (-----)  {:>6.3f} (-----)  {:>6.3f} (-----)    {:<2}\n'
                fmt_2 = '{:<11}  {:>1}     {:>6.3f} (-----)  {:>6.3f} (-----)  {:>6.3f} (-----)         (-----)  {:>6.3f} (-----)    {:<2}\n'
                
                for r in range(1, nr - 1):
                    if id_table[r - 1] >= n1: # columns 3,5,7,11,13 in matlab = 2,4,6,10,12 in python
                        fid.write(fmt_2.format(str(SAtable[r][0]), int(SAtable[r][1]), safe_float(SAtable[r][2]), safe_float(SAtable[r][4]), safe_float(SAtable[r][6]), safe_float(SAtable[r][10]), safe_int(SAtable[r][12])))
                    else: # columns 3,5,7,9,11,13 in matlab = 2,4,6,8,10,12 in python 
                        fid.write(fmt_1.format(str(SAtable[r][0]), int(SAtable[r][1]), safe_float(SAtable[r][2]), safe_float(SAtable[r][4]), safe_float(SAtable[r][6]), safe_float(SAtable[r][8]), safe_float(SAtable[r][10]), safe_int(SAtable[r][12])))
            else:
                fmt_1 = '{:<11}  {:>1}     {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})    {:<2}\n'
                fmt_2 = '{:<11}  {:>1}     {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})                  {:>6.3f} (\261{:>.2f})    {:<2}\n'

                for r in range(1, nr - 1):
                    if id_table[r - 1] >= n1: # columns 3:8,11:13 in matlab = 2:7,10,11,12 in python
                        fid.write(fmt_2.format(str(SAtable[r][0]), int(SAtable[r][1]), safe_float(SAtable[r][2]), abs(safe_float(SAtable[r][3])), safe_float(SAtable[r][4]), abs(safe_float(SAtable[r][5])), safe_float(SAtable[r][6]), abs(safe_float(SAtable[r][7])), safe_float(SAtable[r][10]), abs(safe_float(SAtable[r][11])), safe_int(SAtable[r][12])))
                    else: # columns 3:13 in matlab = 2:12 in python
                        fid.write(fmt_1.format(str(SAtable[r][0]), int(SAtable[r][1]), safe_float(SAtable[r][2]), abs(safe_float(SAtable[r][3])), safe_float(SAtable[r][4]), abs(safe_float(SAtable[r][5])), safe_float(SAtable[r][6]), abs(safe_float(SAtable[r][7])), safe_float(SAtable[r][8]), abs(safe_float(SAtable[r][9])), safe_float(SAtable[r][10]), abs(safe_float(SAtable[r][11])), safe_int(SAtable[r][12])))

            fid.write('------------------------------------------------------------------------------------------------------------ \n')
            
            if K == 1: # columns 3,5,7,11 in matlab = 2,4,6,10 in python
                fmt = '{:<11}  {:>1}     {:>6.3f} (-----)  {:>6.3f} (-----)  {:>6.3f} (-----)         (-----)  {:>6.3f} (-----)\n'
                fid.write(fmt.format(str(SAtable[nr-1][0]), '', safe_float(SAtable[nr-1][2]), safe_float(SAtable[nr-1][4]), safe_float(SAtable[nr-1][6]), safe_float(SAtable[nr-1][10])))
            else: # columns 3:8,11:12 in matlab = 2:7,10,11 in python
                fmt = '{:<11}  {:>1}     {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})  {:>6.3f} (\261{:>.2f})                  {:>6.3f} (\261{:>.2f}) \n'
                fid.write(fmt.format(str(SAtable[nr - 1][0]), '', safe_float(SAtable[nr-1][2]), abs(safe_float(SAtable[nr-1][3])), safe_float(SAtable[nr-1][4]), abs(safe_float(SAtable[nr-1][5])), safe_float(SAtable[nr-1][6]), safe_float(SAtable[nr-1][7]), safe_float(SAtable[nr-1][10]), safe_float(SAtable[nr-1][11])))

            fid.write('============================================================================================================ \n')
            fid.write(' S^a: Structural sensitivity index of individual terms\n')
            fid.write(' S^b: Correlative sensitivity index of individual terms\n')
            fid.write(' S: Total sensitivity index of individual terms\n')
            fid.write(' ST: Total sensitivity index\n')
            if K == 1:
                fid.write(' (--): Cannot compute confidence intervals of listed statistics with K = 1\n')
            else:
                fid.write(f' (\261): {int(100*(1-alfa))}% confidence intervals derived from bootstrapping\n')
                fid.write(' V(i)/V_Y: Relative contribution of each term to model output variance ( = var(Y))\n')
            if K == 1:
                fid.write(' #select: 0 (if term is insignificant) and 1 (significant term)\n')
            else:
                fid.write(' #select: Number of bootstrap trials that identifies respective term as significant\n')

            fid.write('\n')

            fid.write('\n')

    if sys.platform == "win32" or sys.platform == "darwin":
        os.system("start HDMR_settings.txt")
        os.system("start HDMR_results.txt")

    return SA, SA_sig, Fx


def safe_int(value):
    return int(value) if value else 0  # or return a default value if empty


def safe_float(value):
    return float(value) if value else 0.0  # or return a default value if empty


def construct_tables(SA, nterms, p0, RMSE, alfa, f_ret, K, C2, C3, n1, n2, n3, n):
    # This function returns a Table with a) sensitivity, and b) emulator results

    # ----------------------------------------------------------------------- #
    #                 FIRST POSTPROCESS SENSITIVITY ESTIMATES                 #
    # ----------------------------------------------------------------------- #

    # Compute average sensitivity values
    S_m = np.mean(SA['S'], axis=1)
    Sa_m = np.mean(SA['Sa'], axis=1)
    Sb_m = np.mean(SA['Sb'], axis=1)
    
    # Now calculate sum of each statistic
    S_sum = np.sum(SA['S'], axis=0)
    Sa_sum = np.sum(SA['Sa'], axis=0)
    Sb_sum = np.sum(SA['Sb'], axis=0)
    
    # Now calculate associated std's
    Z = lambda p: -np.sqrt(2) * erfcinv(p * 2)
    m = Z(1 - alfa/2)  # Multiplier, alfa is significance level   
    
    # Compute output statistics for Y (variance)
    V_em_div_V_y = SA['V_em'] / SA['V_y']
    V_em_rat = np.mean(V_em_div_V_y, axis=1)
    V_em_div_V_Y_sum = np.sum(V_em_div_V_y, axis=0)
    
    # Compute standard deviation of bootstrap results
    if K > 1:
        S_range = m * np.std(SA['S'], axis=1, ddof=0)
        Sa_range = m * np.std(SA['Sa'], axis=1, ddof=0)
        Sb_range = m * np.std(SA['Sb'], axis=1, ddof=0)
        V_em_range = m * np.std(V_em_div_V_y, axis=1, ddof=0)
        S_sum_range = m * np.std(S_sum)
        Sa_sum_range = m * np.std(Sa_sum)
        Sb_sum_range = m * np.std(Sb_sum)
        V_em_div_V_Y_sum_range = m * np.std(V_em_div_V_Y_sum)
    else:
        S_range, Sa_range, Sb_range, V_em_range = [np.nan] * n, [np.nan] * n, [np.nan] * n, [np.nan] * n

    ST, ST_range = [np.nan] * n1, [np.nan] * n1

    # Now compute the total sensitivity of each parameter/coefficient
    for r in range(n1):
        ij = n1 + np.where(np.sum(C2 == r, axis=1) == 1)[0]
        ijk = n1 + n2 + np.where(np.sum(C3 == r, axis=1) == 1)[0]
        
        # use all bootstrap trials to determine total sensitivity + range!
        TS = np.sum(SA['S'][np.r_[r, ij, ijk], :K], axis=0)
        ST[r] = np.mean(TS)
        if K > 1:
            ST_range[r] = m * np.std(TS)

    # ----------------------------------------------------------------------- #
    #               NOW CONSTRUCT TABULATED TABLE WITH RESULTS                #
    # ----------------------------------------------------------------------- #

    # how many rows of this table
    nr = n1 + n2 + n3 + 1
    
    # initialize row_names of Table and f_i ( = order)
    row_names = [''] * nr
    f_ord = np.full(nr, np.nan)
    
    # now create row_names + f_i
    for r in range(n1):
        f_ord[r] = 1
        row_names[r] = f'x{r+1}'

    for i in range(n2):
        r = i + n1
        f_ord[r] = 2
        row_names[r] = f'x{C2[i, 0]+1}/x{C2[i, 1]+1}'

    for i in range(n3):
        r = i + n1 + n2
        f_ord[r] = 3
        row_names[r] = f'x{C3[i, 0]+1}/x{C3[i, 1]+1}/x{C3[i, 2]+1}'
   
    # add as last row name the sum of the previous rows
    row_names[nr-1] = 'sum'
    
    # now create column names
    col_names = ['term', 'order', 'S_a', 'std.', 'S_b', 'std.', 'S', 'std.', 'S_T', 'std.', 'V(i)/V(Y)', 'std.', '#select']
    nc = len(col_names)

    # Reinitialize SA to become a list of sensitivity analysis results
    SA_table = [[''] * nc for _ in range(nr+1)]
    SA_table[0][:] = col_names
    
    # first column stores row_names
    for i in range(nr):
        SA_table[i+1][0] = row_names[i]

    # Fill columns 2 - 8 of table -> [ f_ord, S_a, std., S_b, std., S, std. ]
    j = [0, 1, 2, 3, 4, 5]
    for r in range(1, nr):
        SA_table[r][1] = int(f_ord[r-1])
        data = np.array([Sa_m[r-1], Sa_range[r-1], Sb_m[r-1], Sb_range[r-1], S_m[r-1], S_range[r-1]])
        for idx in range(len(j)):
            SA_table[r][2 + idx] = f'{data[j[idx]]:4.3f}'

    # Fill columns 9 - 10 of table -> [ ST, std. ]
    j = [0, 1]
    for r in range(1, n1+1):
        data = np.array([ST[r-1], ST_range[r-1]])
        for idx in range(len(j)):
            SA_table[r][8 + idx] = f'{data[j[idx]]:4.3f}'

    # Fill columns 11 - 12 of table -> [ V(i)/V_Y, std.]
    j = [6, 7]
    for r in range(1, nr):
        data = np.array([V_em_rat[r-1], V_em_range[r-1]])
        for idx in range(len(j)):
            SA_table[r][10 + idx] = f'{data[j[idx]-6]:4.3f}'

    # Fill column 13 of table -> [ FX.select ]
    for r in range(1, nr):
        SA_table[r][12] = int(f_ret[r-1])

    # Fill the last row of table (sum)
    if K > 1:
        in_T = [
            np.sum(Sa_m), Sa_sum_range, np.sum(Sb_m), Sb_sum_range, np.sum(S_m), S_sum_range,
            np.nan, np.nan, np.sum(V_em_rat), V_em_div_V_Y_sum_range, np.nan
        ]
    elif K == 1:
        in_T = [
            np.sum(Sa_m), np.nan, np.sum(Sb_m), np.nan, np.sum(S_m), np.nan, np.nan, np.nan,
            np.sum(V_em_rat), np.nan, np.nan
        ]
    
    j = 2 + np.where(~np.isnan(in_T))[0]
    for idx in range(len(j)):
        SA_table[nr][j[idx]] = f'{in_T[j[idx]-2]:4.3f}'

    # ----------------------------------------------------------------------- #
    #               NOW CONSTRUCT TABULATED RESULTS EMULATORS                 #
    # ----------------------------------------------------------------------- #

    row_names = [str(r+1) for r in range(K)]
    col_names = ['emulator', '# terms', '# coefs.', 'RMSE']
    nc = len(col_names)

    Fx = [[''] * nc for _ in range(K+1)]
    Fx[0][:] = col_names
    j = [0, 1]
    for k in range(1, K+1):
        Fx[k][0] = int(k)
        data = np.array([nterms[k-1], p0[k-1]])
        for idx in range(len(j)):
            Fx[k][idx+1] = int(data[j[idx]])

        Fx[k][nc-1] = f'{RMSE[k-1]:4.3f}'

    return SA_table, Fx


def HDMR_plot(SA_sig, Fx, y, y_e, select, p0, iter, id, R, K, refit, maxorder):

    # Output writing to the screen
    print('\n')
    print('HDMR PLOTTING: PLEASE WAIT ...')

    n_program = 'HDMR'
    file_name = f'{n_program}_figures.pdf'
    table1_name = f'Table_1_{n_program}.pdf'
    table2_name = f'Table_2_{n_program}.pdf'
    id_pdf = range(1, table1_name.index('pdf') - 1)

    monitor = get_monitors()[0]
    screen_width = monitor.width
    screen_height = monitor.height
    x_mult = screen_width / 1920
    y_mult = screen_height / 1080
    t_mult = min(x_mult, y_mult)
    fontsize_xylabel = 16 * t_mult
    fontsize_axis = 16 * t_mult
    fontsize_legend = 14 * t_mult
    fontsize_text = 14 * t_mult
    fontsize_table = 16 * t_mult
    fontsize_titlepage = 20 * t_mult

# JAV Added PdfPages to work with Append 
    with PdfPages(file_name) as pdf:

        ### Plot Empty Figure for PDF
        plt.figure(figsize=(12, 6))
        plt.plot([], [], 'ro')  # Empty plot
        plt.axis([0, 1, 0, 1])
        plt.gca().set_facecolor('w')
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])        
        plt.gca().set_xticks([])
        plt.gca().set_yticks([])
        plt.text(0.3 * x_mult, 0.6 * y_mult, r'${\rm Visual \; results \; of \; HDMR \; toolbox}$', fontsize=fontsize_titlepage) #, ha='center', va='center')
        plt.text(0.29 * x_mult, 0.5 * y_mult, r'${\rm Tables \; are \; printed \; to \; this \; PDF \; file}$', fontsize=fontsize_titlepage) #, ha='center', va='center') #, fontweight='bold')
        ax = plt.gca()  # Get current axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        pdf.savefig()    
        plt.show()

        # Plot results of emulator of K trials
        y_plot = np.linspace(np.min(y), np.max(y))
        row = 2
        col = 4

        for k in range(K):
            if k % (row * col) == 0:
                plt.figure(figsize=(15, 10))
                c = 1
                ii = 1
            else:
                ii = ii+1

            if c <= col:
                r = 1
                plot_c = c
            else:
                r = 2
                plot_c = c - col

            id_R = id[0:R, k]
            ax1 = plt.subplot(row, col, ii)
            ax1.plot(y[id_R], y_e[0:R, k], 'rs', markerfacecolor='r', markeredgecolor='r')
            ax1.plot(y_plot, y_plot, '-', color=[0.5, 0.5, 0.5], linewidth=2)
            ax1.set_xlabel(r"$y$", fontsize=fontsize_xylabel)
            if k % col == 0:
                ax1.set_ylabel(r"$y = f(\mathbf{x})$", fontsize=fontsize_xylabel)
            ax1.legend(['Training', '1:1 Line'], loc='lower right', fontsize=fontsize_legend, frameon=False)
            ax1.tick_params(axis='both', labelsize=fontsize_axis)

            RMSE_train = np.sqrt(1 / R * np.sum((y[id_R] - y_e[0:R, k]) ** 2))
            # RMSE_1 = f'{RMSE_train:.4f}'
            ax1.set_title(f'Bootstrap trial: {k + 1}', fontweight='bold', fontsize=fontsize_text)

            ax1.text(0.05, 0.92, f'$\# terms.: {int(np.sum(select[:, k]))}$', transform=ax1.transAxes, fontsize=fontsize_text)
            ax1.text(0.05, 0.85, f'$\# coef.: {int(p0[k])}$', transform=ax1.transAxes, fontsize=fontsize_text)
            ax1.text(0.05, 0.78, f'$\# iter.: {int(np.sum(iter[0:4, k]))}$', transform=ax1.transAxes, fontsize=fontsize_text)
            ax1.text(0.05, 0.70, r"$\text{RMSE}_\text{T}:$"f'{RMSE_train:.3f}', transform=ax1.transAxes, fontsize=fontsize_text)
            c += 1
            if ii == row*col or k == K-1:
                pdf.savefig()
                plt.close()

        plt.close()

        # Plot number of required iterations
        ct = 0
        if K > 1:
            plt.figure(figsize=(10, 5))
#            plt.plot(range(1, K + 1), iter[0, 0:K], 'rs', linewidth=2)
#            plt.plot(range(1, K + 1), iter[0, 0:K], 'rs', markerfacecolor='r', markeredgecolor='r')
            plt.plot(range(1, K + 1), iter[0, 0:K], linestyle='-', linewidth=2, marker='s', color='red', markersize=8, markerfacecolor='red', markeredgewidth=2)
            ct = ct + 1
            if maxorder > 1:
  #              plt.plot(range(1, K + 1), iter[1, 0:K], 'b-', linewidth=2)
  #              plt.plot(range(1, K + 1), iter[1, 0:K], 'bo', markerfacecolor='b', markeredgecolor='b')
                plt.plot(range(1, K + 1), iter[1, 0:K], linestyle='-', linewidth=2, marker='o', color='blue', markersize=8, markerfacecolor='blue', markeredgewidth=2)
                ct = ct + 1
            if maxorder == 3:
#                plt.plot(range(1, K + 1), iter[2, 0:K], 'g-', linewidth=2)
#                plt.plot(range(1, K + 1), iter[2, 0:K], 'gd', markerfacecolor='g', markeredgecolor='g')
                plt.plot(range(1, K + 1), iter[2, 0:K], linestyle='-', linewidth=2, marker='d', color='green', markersize=8, markerfacecolor='green', markeredgewidth=2)
                ct = ct + 1
            if refit == 1:
#                plt.plot(range(1, K + 1), iter[3, 0:K], 'k-', linewidth=2)
#                plt.plot(range(1, K + 1), iter[3, 0:K], 'k^', markerfacecolor='k', markeredgecolor='k')
                plt.plot(range(1, K + 1), iter[3, 0:K], linestyle='-', linewidth=2, marker='^', color='black', markersize=8, markerfacecolor='black', markeredgewidth=2)
                ct = ct + 1
            plt.xlabel('Emulator (bootstrap trial)', fontsize=fontsize_xylabel)
            plt.ylabel('Number of iterations', fontsize=fontsize_xylabel)
            legend_name = ['first order', 'second order', 'third order', 'refitting']
            color_plotting = ['red', 'blue', 'green', 'black']
            # Create custom legend handles with colored lines
            legend_handles = [Line2D([0], [0], color=color_plotting[i], lw=4, label=f"{legend_name[i]}") for i in range(ct)]
            legend = plt.legend(handles=legend_handles, fontsize=fontsize_legend, loc='best')
            # Match text color to line color
            for z, text in enumerate(legend.get_texts()):
                text.set_color(color_plotting[z])  
        
#            plt.legend(['first order', 'second order', 'third order', 'refitting'], fontsize=fontsize_legend)

            plt.tight_layout()
            pdf.savefig()
            plt.close()

        # Plot Fx to table in figure
        Fx = np.array(Fx)
        Fx_df = pd.DataFrame(Fx[1:K + 1, 0:4], columns=[Fx[0, i] for i in range(0,4)])
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.axis('tight')
        ax.axis('off')
        ax.table(cellText=Fx_df.values, colLabels=Fx_df.columns, loc='center', cellLoc='center', colLoc='center')
        plt.title('Emulator performance training data set', fontsize=fontsize_titlepage)
        # plt.savefig(table1_name, format='pdf')
        pdf.savefig()
        plt.close()

        # Plot SA_sig to table in figure
        SA_sig = np.array(SA_sig)
        # SA_sig_df = pd.DataFrame(SA_sig, columns=[SA_sig[0, i] for i in range(0, SA_sig.shape[1])])    used in HDMR_EXT
        SA_sig_df = pd.DataFrame(SA_sig[1:, 1:], columns=[SA_sig[0, i] for i in range(1, SA_sig.shape[1])])
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.axis('tight')
        ax.axis('off')
        ax.table(cellText=SA_sig_df.values, colLabels=SA_sig_df.columns, loc='center', cellLoc='center', colLoc='center')
        plt.title('Variance-based decomposition and sensitivity coefficients', fontsize=fontsize_titlepage)
        # plt.savefig(table2_name, format='pdf')
        pdf.savefig()
        plt.close()

    # Open the final PDFs
    os.startfile(file_name)
    #os.startfile(table1_name)
    #os.startfile(table2_name)

# Latin Hypercube Sampling function
def LH_sampling(mn, mx, N):
    """
    Latin Hypercube Sampling.
    
    Args:
        mn: Lower bound vector
        mx: Upper bound vector
        N: Number of samples to generate

    Returns:
        N x d matrix of Latin Hypercube samples
    """
    d = len(mn)  # Number of parameters
    rng = np.array(mx) - np.array(mn)  # 1 x d vector with parameter ranges
    y = np.random.rand(N, d)  # N x d matrix with uniform random labels
    # id_matrix = np.argsort(np.random.rand(N, d), axis=0)  # Random sort (1:N without replacement)
    # really important change below so that X stays in bound! as list is from 0 - N-1 rather than 1 to N
    id_matrix = 1 + np.argsort(np.random.rand(N, d), axis=0)  # Random sort (1:N without replacement)
    M = (id_matrix - y) / N  # Multiplier matrix (y introduces randomness)
    R = np.add(np.multiply(M, rng), mn)  # N x d matrix of stratified LH samples
    
    return R
