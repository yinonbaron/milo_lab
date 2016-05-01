# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 09:23:15 2016

@author: yinonbaron
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def f(y, t, params):
    a, R, A = y                         # unpack current values of y
    K_a, K_r, v, gamma, mu = params     # unpack parameters
    derivs = [v(a)*A-gamma(a)*R,        # da/dt
             K_r(a)-mu(a,R)*R,            # dR/dt
             K_a(a)-mu(a,R)*A             # dA/dt
             ]
    return derivs

# Parameters
N_prot = 3*1e6
aa_prot = 300.0
g0 = 25.0
K_g = 1.0

v0 = 20.0
K_v = 4.0

k_A0 = 10.0
K_A = 1.0

k_R0 = 10.0
K_R = 1.0
K_a = lambda x: k_A0/(1+(x/K_A)**2)
K_r = lambda x: (k_R0*x)/(K_R+x) # MM kinetics for the activation of the ribosomal rRNA as a function of amino-acid concentration
v = lambda x: v0/(1+(x/K_v)**2)
gamma = lambda x: g0*((x/K_g)**2/(1+(x/K_g)**2))
mu = lambda x, R: np.log(2)/((N_prot*aa_prot)/(gamma(x)*R)) # ln2/ time it takes in seconds to translate the total amount of aa in a cell

# Initial values
a0 = 0.0     # initial angular displacement
R0 = 0.0     # initial angular velocity
A0 = 0.0     # initial angular displacement
# Bundle parameters for ODE solver
params = [K_a, K_r, v, gamma, mu]

# Bundle initial conditions for ODE solver
y0 = [a0, R0, A0]

# Make time array for solution
tStop = 200.
tInc = 0.05
t = np.arange(0., tStop, tInc)

# Call the ODE solver
psoln = odeint(f, y0, t, args=(params,))