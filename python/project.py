import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, exp, sin, cos
import time
from scipy.optimize import minimize
from flight_time import obj_func, obj_func_print, scale, unscale, setPayload


if __name__ == '__main__':

    from pyoptwrapper import optimize
    from pyoptsparse import NSGA2, SNOPT, NLPQLP, ALPSO, SLSQP

    eng_displace = 80.0 # cc
    fuel_rate = .02 # kg/s
    omega_r = 1500 # rpm
    R_rotor = 0.705 # m
    fuel_cap = 1.0 # kg
    theta_hover = 20 # deg
    x0 = scale([eng_displace, fuel_rate, omega_r, R_rotor, fuel_cap, theta_hover])
    # x0 = [ 4.3441622,   0.16620041,  0.10885692,  2.78460696,  6.76453586, 0.85782237]
    x0 = [ 4.61280557,  0.17376677,  0.08200489,  3.4691025,   7.32432873, 0.86800068]

    lb = [0, 0, 0, 0, 0, 0]
    ub = scale([1200.0, 1.0, 7500.0, 3.0, 12.0, 30.0])
    setPayload(5.0);

#    print obj_func_print(x0)

    optimizer = SNOPT()
    # optimizer = NSGA2()
    # optimizer.setOption('maxGen', 200)

    xopt, fopt, info = optimize(obj_func, x0, lb, ub, optimizer)
    print 'SNOPT:', xopt, fopt, info

    out = obj_func_print(xopt)

    steps = 30
    payloads = np.linspace(0.0, 8.0, num=steps)
    flight_times = np.zeros(steps)
    engine_size = np.zeros(steps)
    for i in range(steps):
        setPayload(payloads[i]);
        x0 = xopt
        xopt, fopt, info = optimize(obj_func, x0, lb, ub, optimizer)
        flight_times[i] = -fopt/60
        engine_size[i] = xopt[0]*80.0

    plt.figure()
    plt.xlabel('Useful Payload (kg)')  # labels for axes
    plt.ylabel('Flight Time (minutes)')
    plt.title('Pareto Front')
    plt.plot(payloads, flight_times, '-k', label='Data')
    # plt.legend(loc='upper left')
    # plt.show()

    plt.figure(2)
    plt.xlabel('Useful Payload (kg)')  # labels for axes
    plt.ylabel('Engine Size (cc)')
    plt.title('Engine Size with Payload')
    plt.plot(payloads, engine_size, '-k', label='Data')
    # plt.legend(loc='upper left')
    plt.show()
