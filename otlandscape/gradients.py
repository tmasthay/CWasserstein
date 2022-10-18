import numpy as np
import copy

def numerical_deriv(x, f, tol):
    n = len(x)
    y = np.zeros((n,))
    for i in range(n):
        perturb = copy.copy(x)
        perturb[i] = x[i] + tol
        left = f(perturb)
        right = f(x)
        y[i] = (1.0 / tol) * (left - right)
        print('GRADIENT NUMERATOR: %s ... ... ... GRADIENT: %s'%((left - right), y))
    return y

def numerical_hessian(x, f, tol):
    n = len(x)
    y = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if( i != j ):
                pll = copy.copy(x)
                plr = copy.copy(x)
                prl = copy.copy(x)
                prr = copy.copy(x)

                pll[i] = x[i] + tol
                pll[j] = x[j] + tol
 
                plr[i] = x[i] + tol
                plr[j] = x[i] - tol
 
                prl[i] = x[i] - tol
                prl[j] = x[j] + tol
 
                prr[i] = x[i] - tol
                prr[j] = x[j] - tol
 
                y[i][j] = (0.25 / tol**2) * (f(pll) - f(plr) - f(prl) + f(prr))

            else:
                forward = copy.copy(x)
                backward = copy.copy(x)
                forward[i] = x[i] + tol
                backward[i] = x[i] - tol
                y[i][j] = (0.5 / tol**2) * (f(forward) - 2 * f(x) + f(backward))
    return y
        
     

