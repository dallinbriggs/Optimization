import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, exp, sin, cos
import time
from scipy.optimize import minimize
from flight_time import obj_func, obj_func_print, scale, unscale


if __name__ == '__main__':

    from pyoptwrapper import optimize
    from pyoptsparse import NSGA2, SNOPT, NLPQLP, ALPSO, SLSQP

    eng_displace = 80.0 # cc
    fuel_rate = .02 # kg/s
    omega_r = 1500 # rpm
    R_rotor = 0.705 # m
    c_rotor = 0.05 # m
    fuel_cap = 1.0 # kg
    theta_hover = 20 # deg
    x0 = scale([eng_displace, fuel_rate, omega_r, R_rotor, c_rotor, fuel_cap, theta_hover])

    lb = [0, 0, 0, 0, 0, 0, 0]
    ub = scale([1200.0, 1.0, 7500.0, 3.0, 0.1, 12.0, 30.0])

#    print obj_func_print(x0)

    optimizer = SNOPT()
#    optimizer = NSGA2()
#    optimizer.setOption('maxGen', 200)

    xopt, fopt, info = optimize(obj_func, x0, lb, ub, optimizer)
    print 'SNOPT:', fopt, info

    out = obj_func_print(xopt)



