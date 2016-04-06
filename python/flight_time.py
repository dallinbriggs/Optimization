#!/usr/bin/env python
import numpy as np
from math import sqrt, exp, sin, cos, pi

prt = False

def flight_time(eng_displace, fuel_rate, omega, R_rotor, c_rotor, f_cap, theta_hover, n_rotors, payload, frame, rho_air, E_rhof, n_engine):

#    --- design variables ---
#    eng_displace    engine displacement (cubic centemeter)
#    fuel_rate       rate that fuel is consumed (kg/s)
#    omega           angular velocity of the rotor (rad/s)
#    R_rotor         rotor radius (meters)
#    c_rotor         cord of the rotor (meters)
#    f_cap           fuel capacity (liters)
#    theta_hover     pitch angel of the blade in hover (rad)
#    --- parameters ---
#    n_rotors        number of 2 blade rotors  !!!! this is a discrete variable that can't be dealt with yet but could be later
#    payload         mass of the payload (kg)
#    frame           mass of the frame (kg) !!!! this could be calulated from all up weight and radius of the rotor
#    rho_air         air density (kg/m^3)
#    E_rhof          energy density of fuel (Joules/kg)
#    n_engine        efficiency of the engine !!!! this is a param for now but could be calculated from an engine model

    W_motor = eng_displace*0.02145 + 0.25412   #mass of the motor (Kg) !!!! this is an emperical best fit linear relationship
    P_max = eng_displace*76.0 + 502.0 #max power (W) as a function of displacment based on emperical best fit linear relationship
    P_max = -0.0454790393518*eng_displace**2 + 79.5042893755*eng_displace #max power (W) as a function of displacment, best fit linear quadradic
    W_f = f_cap*.75                            #Fuel mass (kg)
    AUW = W_motor + W_f + frame + payload      #All up mass in kg (about 52 lbs)

#    --- first find a induced air velocity ---  !!! this is an assumption that we will need to improve later
    T = AUW*9.8/n_rotors                       #Thrust required
    V_i = sqrt(T/(2*pi*rho_air*R_rotor**2))     #Induced air velocity

#   --- next integrate over the blade to produce T thrust and Q torque ---
#    http://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node86.html

    alpha0 = 0.0
    CL_0 = 0.28
    CL_alpha = 3.45
    M = 0.55
    CD_p = 0.0437

    T = 0.0
    Q = 0.0

    for i in range(10):  # 10 stations...
        dr = R_rotor/10
        Vn = omega*dr*(i+1)
        Ve = sqrt(Vn**2 + V_i**2)
        alpha = theta_hover - np.arctan(V_i/Vn)
        phi = np.arctan(V_i/Vn)

        sigma_a = (1 + exp(-(M*(alpha - alpha0))) + exp((M*(alpha + alpha0))))/((1 + exp(-(M*(alpha - alpha0))))*(1 + exp((M*(alpha + alpha0))))) # to model stall
        CL_a = (1 - sigma_a)*(CL_0 + CL_alpha*alpha) + sigma_a*(2*np.sign(alpha)*pow(sin(alpha),2.0)*cos(alpha))
        AR = 2*R_rotor/c_rotor
        CD_a = CD_p + ((pow((CL_0 + CL_alpha*(alpha)),2.0))/(3.14159*0.9*AR))

#        differential drag and lift
        dL = 0.5*rho_air*(Ve**2)*c_rotor*CL_a*dr
        dD = 0.5*rho_air*(Ve**2)*c_rotor*CD_a*dr

#        differential thrust and torque
        dT = dL*cos(phi + alpha) - dD*sin(phi + alpha)
        dQ = dr*(i+1)*(dL*sin(phi + alpha) + dD*cos(phi + alpha))

#        integrate over the length of the rotor
        T += dT
        Q += dQ
    T = n_rotors*2*T
    Q = n_rotors*2*Q

#    T, Q = nasa(2.0, R_rotor, omega, n_rotors, theta_hover, AUW, rho_air)

    PowerRequired = 8*Q*omega
    PowerProduced = fuel_rate*E_rhof*n_engine
    c = np.array([PowerProduced - PowerRequired, T - AUW*9.8 , P_max - PowerProduced])
    global prt
    if prt == True:
        print 'PP', PowerProduced
        print 'PR', PowerRequired
        print 'T', T
        print 'W', AUW*9.8
        print 'Pmax',P_max
    FT = W_f/fuel_rate
    dft_dx = np.array([0, -W_f/(fuel_rate**2), 0, 0, 0, .75, 0])
#    dft_dx = np.array([-W_f/(fuel_rate**2), 0, 0])
    return FT, c, dft_dx, []

def obj_func(x):

    eng_displace, fuel_rate, omega_r, R_rotor, c_rotor, fuel_cap, theta_hover = unscale(x)

    omega_r = omega_r*2*3.14159/60 # to rad/s
    theta_hover = theta_hover*3.14159/180 # to rad

    n_rotors = 6
    payload = 3.0
    frame = 2.0
    rho_air = 1.225
    E_rhof = 44.4e6
    n_engine = 0.15

    ft, c, dft_dx, dc_dx = flight_time(eng_displace, fuel_rate, omega_r, R_rotor, c_rotor, fuel_cap, theta_hover, n_rotors, payload, frame, rho_air, E_rhof, n_engine)
    return -ft, -c

def obj_func_print(x):

    global prt
    prt = True
    ft, c = obj_func(x)
    prt = False

    eng_displace, fuel_rate, omega_r, R_rotor, c_rotor, fuel_cap, theta_hover = unscale(x)

    print 'displacement (cc)',eng_displace
    print 'fuel rate (kg/s)',fuel_rate
    print 'angular velocity of the rotor (rpm)',omega_r
    print 'Rotor radius (m)',R_rotor
    print 'rotor chord (m)',c_rotor
    print 'fuel capacity (kg)',fuel_cap
    print 'theta_hover (deg)',theta_hover
    print 'flight time (min)',-ft/60.0
    return ft, c

def scale(x):
    eng_displace, fuel_rate, omega_r, R_rotor, c_rotor, fuel_cap, theta_hover = x

    eng_displace = eng_displace/80.0
    fuel_rate = fuel_rate/.02
    omega_r = omega_r/1500.0
    R_rotor = R_rotor/0.705
    c_rotor = c_rotor/0.05
    fuel_cap = fuel_cap
    theta_hover = theta_hover/20

    return [eng_displace, fuel_rate, omega_r, R_rotor, c_rotor, fuel_cap, theta_hover]

def unscale(x):
    eng_displace, fuel_rate, omega_r, R_rotor, c_rotor, fuel_cap, theta_hover = x

    eng_displace = eng_displace*80.0
    fuel_rate = fuel_rate*.02
    omega_r = omega_r*1500.0
    R_rotor = R_rotor*0.705
    c_rotor = c_rotor*0.05
    fuel_cap = fuel_cap
    theta_hover = theta_hover*20

    return [eng_displace, fuel_rate, omega_r, R_rotor, c_rotor, fuel_cap, theta_hover]

'''
def obj_func_penelty(x):

    eng_displace, fuel_rate, omega_r, R_rotor, c_rotor, fuel_cap, theta_hover = x
#    fuel_rate, R_rotor = x

    eng_displace = eng_displace*80.0
    fuel_rate = fuel_rate*.02
    omega_r = omega_r*1500*2*3.14159/60 #rpm to rad/s
    R_rotor = R_rotor*0.705
    c_rotor = c_rotor*0.05
    fuel_cap = fuel_cap*1.0
    theta_hover = theta_hover*6*3.14159/180

    n_rotors = 4
    payload = 5.0
    frame = 2.0
    rho_air = 1.225
    E_rhof = 44.4e6
    n_engine = 0.15

    ft, c, dft_dx, dc_dx = flight_time(eng_displace, fuel_rate, omega_r, R_rotor, c_rotor, fuel_cap, theta_hover, n_rotors, payload, frame, rho_air, E_rhof, n_engine)

    P = 0.0;
    for i in range(len(c)):
        if c[i] > 0:
            P += 100.0/2.0*c[i]**2

    return -ft + P


def nasa(N, R, omega, n, theta, AUW, rho):

#    N           Number of blades in each rotor. A design variable?
#    R           Rotor radius in meters, a design variable.
#    omega       rotor system angular velocity, a design variable.
#    n           Number of rotors, tri, vs quad, etc.
#    theta       Blade commanded pitch angle, in rad.
#    AUW         All up weight (kg)
#    rho         Air density.
    a = pi**2/90                               #Blade lift-curve slope. Typical value for thin airfoils in small angles of attack from wikipedia.
    c = .09*R                                  #Blade chord in meters, probably a function of blade radius. Based on typical rc rotor blade.
    e = .05*R                                  #Flapping hing offset in meters, typically 2 to 5 percent the rotor radius. RC helicopters don't have flapping hinges, but they have lead lag hinges which may act like flapping hinges.
    epsilon = e/R
    V = 0                                      #True airspeed, m/sec. I believe this is the airspeed of the whole aircraft.
    alpha = 0                                  #Hub plane angle of attack, also not sure if the unit is in degrees or radians. Paper says deg. 0 because it is in a hover.
    mu = V*cos(alpha)/(omega*R)                #advance ratio. 0 because it is hovering.
    g = 9.8
    v_i = sqrt(AUW*g/(n*2*rho))                #Uniform induced velocity, m/sec. (C_T*omega*R)/(2*(mu^2+lambda^2)^.5) Circular reference with lambda. Current equation is based on induced velocity for one weight.
    lmda = (V*sin(alpha) - v_i)/(omega*R)    #inflow ratio, used in thrust equation.
    theta_0 = theta                            #Blade-root collective pitch measured from hub plane in rad.
    theta_t = 0                                #Total blade twist, the paper says in degrees, but I think it should be in radians. 0 because rc heli rotors have no twist.
    B_lc = 6*pi/180                            #longitudinal cycle pitch measured from the hub plane and in "wind-hub" system in rad (or deg).
    K_l = 0                                    #Pitch-flap coupling ratio, based on most rc helicopter designs.
    b_l = 1*pi/180                             #?? lateral first-harmonic flapping coefficient measured from hub plane in "wind-hub" system, rad.
    a_0 = 1*pi/180                             #Blade coning angle measured from hub plane, rad. Small because rpm is high and blades are short.
    a_l = 1*pi/180                             #?? longitudinal first-harmonic flapping coefficient measuered from hub plane and in the "wind-hub" system, rad.
    adot_0 = 0                                 #derivative of blade coning angle.
    adot_l = 0                                 #derivative of longitudinal first-harmonic flapping coefficient.
    bdot_l = 0                                 #?? derivative of lateral flapping coefficient.
    p = 0                                      #Aircraft roll rate, rad/sec. Used in thrust equation.
    beta_w = 0                                 #Rotor sideslip angle, the angle between x_s and x'_s, in deg. Assume no wind.
    q = 0                                      #Aircraft pitch rate, rad/sec. Used in thrust equation.
    M_beta = (.47*R-.1342)*R**2/3              #Blade mass moment about the flapping hinge, kg-m.
    addot_0 = 0                                #second derivative of blade coning angle measured from hub plane.
    delta = .0141                              #blade mean profile drag coefficient, based on naca 0015 at 5 degrees.
    A_lc = 5*pi/180                            #lateral cyclic pitch measured from hub plane and in "wind-hub" system, in rad (or deg).


    Thrust = N/2*rho*a*c*R*(omega*R)**2*(.5*(1-epsilon**2)*lmda + theta_0*(1/3 + mu**2/2*(1-epsilon)) +  \
             theta_t*(.25+mu**2/2*(1-epsilon**2)) - mu/2*(1-epsilon**2)*(B_lc-K_l*b_l) - a_0*(1/3+mu**2/2*(1-epsilon))*K_l  \
             + a_l*(mu/2*epsilon*(1-epsilon)) - adot_0/omega*(1/3-epsilon/2) + bdot_l/omega*(mu/4*(1-epsilon)**2)  \
             + mu/4*(1-epsilon**2)*(p/omega*cos(beta_w) + q/omega*sin(beta_w)) ) - N*M_beta*addot_0

    Q_coeff = delta/(4*a)*(1 + (1-epsilon**2)*mu**2) - (theta_0 - K_l*a_0)* \
              (lmda/3+(epsilon/3-.25)*adot_0/omega + mu/6*bdot_l/omega \
              - mu*epsilon/4*(bdot_l/omega-a_l) + mu/6*(p/omega*cos(beta_w) \
              + q/omega*sin(beta_w))) - (A_lc-K_l*a_l)*((.125-epsilon/6)*(adot_l/omega+b_l) \
              - mu/6*a_0+b_l/16*(1-epsilon**2)*mu**2+1/8*(-p/omega*sin(beta_w)+q/omega*cos(beta_w))) \
              - (B_lc-K_l*b_l)*((.125-epsilon/6)*(bdot_l/omega-a_l) \
              + (epsilon/4-1/6)*mu*adot_0/omega+.5*(1-epsilon**2)*(mu*lmda/2+a_l/8*mu**2) \
              + 1/8*(p/omega*cos(beta_w)+q/omega*sin(beta_w))) - theta_t*(lmda/4  \
              + (epsilon/4-1/5)*adot_0/omega+mu/8*(bdot_l/omega)-epsilon*mu/6*(bdot_l/omega-a_l))  \
              - .5*(1-epsilon**2)*(lmda**2+lmda*mu*a_l+2*lmda*epsilon*adot_0/omega  \
              + mu*epsilon*(a_l*adot_0/omega+a_0*(adot_l/omega+b_l))  \
              + mu**2*(a_0**2/2+3/8*a_l**2+1/8*b_l**2)) + mu/3*(a_l*(adot_0/omega) \
              + a_0*(adot_l/omega+b_l))+2/3*lmda*(adot_l/omega)  \
              - (-mu/3*a_0 + (.25-epsilon/3)*(adot_l/omega+b_l))* \
              (-p/omega*sin(beta_w)+q/omega*cos(beta_w))  \
              - (.25-epsilon/3)*(bdot_l/omega-a_l)*(p/omega*cos(beta_w) \
              + q/omega*sin(beta_w)) - .125*(-p/omega*sin(beta_w)+q/omega*cos(beta_w))**2  \
              - .125*(p/omega*cos(beta_w)+q/omega*sin(beta_w))**2 -  \
              (.25-2/3*epsilon+epsilon**2/2)*((adot_0/omega)**2  \
              + .5*((adot_l/omega+b_l)**2+(bdot_l/omega-a_l)**2))

    Torque = N/2*rho*a*c*R**2*(omega*R)**2*(Q_coeff)

    return Thrust, Torque
'''
