import numpy as np
from scipy.optimize import minimize

A = np.array([
    [-0.13086780552598704, -0.9913998272527651], 
    [0.72730025766706086, -0.6863194119339967]
  ])

#Epsilon defined in paper
e = 0.1

def cmFunction(x):
    """Function to minimize"""
    n0 = np.array([x[0], x[1], 1])
    n1 = np.array([x[2], x[3], 1])
    
    t10 = np.array([A[1][0], A[1][1], x[5]])
    t01 = np.array([A[0][0], A[0][1], x[4]])
    
    return ((np.power(np.linalg.norm(np.cross(t10, n0)),2)) + (np.power(np.linalg.norm(np.cross(t01, n1)),2)) + np.power(x[4],2) + np.power(x[5],2));
   
cons = ({'type': 'ineq', 'fun': lambda x:  (x[0]*x[2] + x[1]*x[3] + 1.0) + e},
        {'type': 'ineq', 'fun': lambda x: -(x[0]*x[2] + x[1]*x[3] + 1.0) + e},
        {'type': 'ineq', 'fun': lambda x: (A[0][0]*A[1][0] + A[0][1]*A[1][1] + x[4]*x[5]) + e},
        {'type': 'ineq', 'fun': lambda x: -(A[0][0]*A[1][0] + A[0][1]*A[1][1] + x[4]*x[5]) + e},
        {'type': 'eq', 'fun': lambda x:  (A[0][0]*x[0] + A[0][1]*x[1] + x[4])},
        {'type': 'eq', 'fun': lambda x:  (A[1][0]*x[2] + A[1][1]*x[3] + x[5])})

x0 = np.array([1,1, 1,1, 1,1]) 
res = minimize(cmFunction, x0, method='SLSQP', constraints=cons, options={'xtol': 1e-8, 'disp': True})
print res