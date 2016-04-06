import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, exp, sin, cos
import time
from scipy.optimize import minimize
from flight_time import obj_func


if __name__ == '__main__':

    from pyoptwrapper import optimize
    from pyoptsparse import NSGA2, SNOPT, NLPQLP, ALPSO, SLSQP

    eng_displace = 80.0
    fuel_rate = .02
    omega_r = 1500*2*3.14159/60 #rpm to rad/s
    R_rotor = 0.705
    c_rotor = 0.05
    fuel_cap = 1.0
    theta_hover = 20*3.14159/180
    #    x0 = [eng_displace, fuel_rate, omega_r, R_rotor, c_rotor, fuel_cap, theta_hover]
    x0 = [1., 1., 1., 1., 1., 1., 1.]
    #    x0 = [fuel_rate, R_rotor, theta_hover]

    lb = [0, 0, 0, 0, 0, 0, 0]
    ub = [1200.0/80.0, 1.0/.02, 7500.0/1500.0, 2.0/0.705, 0.1/0.05, 12.0/1.0, 30.0/20.0]

#    x0 = [3.55613094,   0.24996519,   0.72821006,   1.27615376,   1.53249666,  12.,
#    1.32161584]
#    x0 = [  9.85947913e-04,   1.73583194e-01,   4.12512324e-01 ,  1.84357306e+00,
#    1.75334613e+00,   1.19999999e+01,   6.94509532e-01]
    x0 = [  3.37894914e-05,   2.46734305e-01,   3.92570916e-01,   1.76376878e+00,
    1.24587979e+00,   1.19999995e+01,   1.47021053e+00]

    print obj_func(x0)

#    optimizer = SNOPT()
    optimizer = NSGA2()
    optimizer.setOption('maxGen', 200)

    xopt, fopt, info = optimize(obj_func, x0, lb, ub, optimizer)
    print ub
    print 'SNOPT:', xopt, fopt, info

    print obj_func(xopt)



