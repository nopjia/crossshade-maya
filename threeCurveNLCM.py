import numpy as np
from scipy.optimize import minimize

A = np.array([[0.98107917124316835, -0.19360697237656041],[0.87671060288001323, 0.48101821046376625],[0.98205077997564239, 0.18861671598570601], [0.91552166056331108, -0.40226867767624847]])

#Epsilon defined in paper
e = 0.1

def cmFunction(x):
    """Function to minimize"""
    n0 = np.array([x[0], x[1], 1])
    n1 = np.array([x[2], x[3], 1])
    
    t10 = np.array([A[1][0], A[1][1], x[7]])
    t01 = np.array([A[0][0], A[0][1], x[6]])
    
    s = ((np.power(np.linalg.norm(np.cross(t10, n0)),2)) + (np.power(np.linalg.norm(np.cross(t01, n1)),2)) + np.power(x[6],2) + np.power(x[7],2))
    n2 = np.array([x[4], x[5], 1])
    
    t21 = np.array([A[3][0], A[3][1], x[9]])         
    t12 = np.array([A[2][0], A[2][1], x[8]])
    
    return s + ((np.power(np.linalg.norm(np.cross(t21, n1)),2)) + (np.power(np.linalg.norm(np.cross(t12, n2)),2)) + np.power(x[8],2) + np.power(x[9],2));
   
cons = ({'type': 'ineq', 'fun': lambda x:  (x[0]*x[2] + x[1]*x[3] + 1.0) + e},
        {'type': 'ineq', 'fun': lambda x: -(x[0]*x[2] + x[1]*x[3] + 1.0) + e},
        {'type': 'ineq', 'fun': lambda x: (A[0][0]*A[1][0] + A[0][1]*A[1][1] + x[6]*x[7]) + e},
        {'type': 'ineq', 'fun': lambda x: -(A[0][0]*A[1][0] + A[0][1]*A[1][1] + x[6]*x[7]) + e},       
        {'type': 'ineq', 'fun': lambda x:  (x[2]*x[4] + x[3]*x[5] + 1.0) + e},
        {'type': 'ineq', 'fun': lambda x: -(x[2]*x[4] + x[3]*x[5] + 1.0) + e},
        {'type': 'ineq', 'fun': lambda x: (A[2][0]*A[3][0] + A[2][1]*A[3][1] + x[8]*x[9]) + e},
        {'type': 'ineq', 'fun': lambda x: -(A[2][0]*A[3][0] + A[2][1]*A[3][1] + x[8]*x[9]) + e},      
        {'type': 'eq', 'fun': lambda x:  (A[0][0]*x[0] + A[0][1]*x[1] + x[6])},
        {'type': 'eq', 'fun': lambda x:  (A[1][0]*x[2] + A[1][1]*x[3] + x[7])},
        {'type': 'eq', 'fun': lambda x:  (A[2][0]*x[2] + A[2][1]*x[3] + x[8])},
        {'type': 'eq', 'fun': lambda x:  (A[3][0]*x[4] + A[3][1]*x[5] + x[9])})


x0 = np.array([-1,-1, -1,-1, -1,-1, -1,-1, -1,-1])
res = minimize(cmFunction, x0, method='SLSQP', constraints=cons, options={'xtol': 1e-8, 'disp': True})
print res