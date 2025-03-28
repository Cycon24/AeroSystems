# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 21:41:37 2023

@author: cycon
"""

import numpy as np
from scipy.linalg import solve


def calc_static_props_1(To, Po, mdot, A, gam=1.4, R=287, initial_g_dec=0.5, maxIter=1000, tolerance=1e-6):
    # Define the functions and their derivatives
    def f(T, P, A, To):
        return T + (T**2 / P**2) * A - To
    
    def g(T, P, A, Po, k):
        return P**k + P**(k-2) * T * A - Po**k
    
    def df_T(T, P, A):
        return 1 + (2 * T / P**2) * A
    
    def df_P(T, P, A):
        return - (2 * T / P**3) * A
    
    def dg_T(T, P, A, k):
        return P**(k-2) * A
    
    def dg_P(T, P, A, k):
        return k * P**(k-1) + (k-2) * P**(k-2) * T * A
    
    # Newton-Raphson iteration
    def newton_raphson(T, P, A, To, Po, k, tol=1e-2, max_iter=100000):
        for i in range(max_iter):
            # Evaluate functions and derivatives
            F = np.array([f(T, P, A, To), g(T, P, A, Po, k)])
            J = np.array([[df_T(T, P, A), df_P(T, P, A)],
                          [dg_T(T, P, A, k), dg_P(T, P, A, k)]])
    
            # Solve for the update
            delta = solve(J, -F)
    
            # Update variables
            T += delta[0]
            P += delta[1]
    
            # Check for convergence
            if np.linalg.norm(delta) <= tol:
                # print(f"Converged after {i+1} iterations.")
                return T, P
    
        # print("Did not converge within the maximum number of iterations.")
        return None
    
    # Known inputs
    A = ((gam-1)*R*mdot**2) / (2*gam*A**2)
    To = To
    Po = Po
    k = (gam-1)/gam
    
    # Initial guess
    # initial_guess = np.array([0.9*To, 0.9*Po])
    
    # Call the Newton-Raphson solver
    ig_start = int(initial_g_dec*1000)
    for p in range(ig_start,1000,1):
        i = p/1000
        initial_guess = np.array([i*To, i*Po])
        result = newton_raphson(initial_guess[0], initial_guess[1], A, To, Po, k, tolerance, maxIter)
        if result:
            print('Initial Guess: {:.2%} of Stagnation'.format(p/1000))
            break
    # Print the result
    if result:
        T_sol, P_sol = result
        rho = P_sol/(R*T_sol)
        print("Solution: \nT = {:.3f} K \nP = {:.3f} Pa \nρ = {:.3f} kg/m^3".format(T_sol, P_sol, rho))
        result = [T_sol, P_sol, rho]
    else:
        print('Could not converge on a solution for Static Properties within combustor.')
        result = [None, None, None]
    return result


def calc_static_props_2(To, Po, V, gam=1.4, R=287):
    T = To - (V**2)*(gam - 1)/(2*gam*R)
    P = Po/((To/T)**(gam/(gam-1)))
    rho = P/(R*T)
    print("Solution: \nT = {:.3f} K \nP = {:.3f} Pa \nρ = {:.3f} kg/m^3".format(T, P, rho))
    return [T, P, rho]

def calc_area(To, Po, V, mdot, gam=1.4, R=287):
    [T, P, rho] = calc_static_props_2(To, Po, V, gam, R)
    Area = mdot/(rho*V)
    return Area
    

# Test
if __name__=='__main__':
    To_i = 475
    Po_i = 4.552e5
    mdot = 9
    Area = 0.0389
    [T, P, rho] = calc_static_props_1(To_i, Po_i, mdot, Area, initial_g_dec=0.98, maxIter=10000, tolerance=1e-2)
    
    V = mdot/(rho*Area)
    [T, P, rho] = calc_static_props_2(To_i, Po_i, V)
    print('Area: {:.4f}'.format(calc_area(To_i, Po_i, V, mdot)))
