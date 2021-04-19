import numpy as np
from scipy.optimize import fsolve
import scipy.integrate as integrate

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

"""
Constants
"""
Vtip = 218
R = 11.94/2
Omega = Vtip/R
gamma = 6
cla = 5.7
sigma = 0.075
rho = 1.225
CD = 1
S = 1.5
g = 9.81
m = 4300
Iy = 10000
h = 1

"""
While loop definition
"""
Tstart = 0
Tstop = 120
N_timesteps = 1200
dt = Tstop/N_timesteps

t = 0

"""
Control inputs
"""
theta_0 = 6*np.pi/180
theta_c = 0

theta_0gen = 5*np.pi/180
"""
State variables
"""
u = 46.3
w = 0
q = 0
theta_f = 0
x = 0
y = 0
"""
Initialize result arrays
"""
lambda_i = np.sqrt(m*g/(2*rho*(np.pi*R**2)))/Vtip; # Initial guess
ulst = []
wlst = []
qlst = []
theta_flst = []
tlst = []
wdot = 0
xlst = []
ylst = []
wdes = 0
ydes = 0
ydotdes = 0
ydot = 0
ydotlst = []
#thetaclst = []

"""
Gains
"""

K1 = 0.006
K2 = 0.0001
K3 = 0.001

t_control = 15   #Time after which the pilot becomes active

while t < Tstop:
    
    if t > t_control:
        theta_c = 0.2*theta_f + 0.2*q

    V = np.sqrt(u**2 + w**2)

    if u == 0:
        if w > 0:
            alpha_c = theta_c + np.pi/2
        else:
            alpha_c = theta_c - np.pi/2
    else:
        alpha_c = theta_c - np.arctan2(w,u)
    
    mu = ((V)/(Omega*R))*np.cos(alpha_c)
    lambda_c = (V*np.sin(alpha_c))/(Omega*R)

    a1 = ((8*mu*theta_0)/(3) - 2*mu*(lambda_c + lambda_i) - (16*q)/(gamma*Omega))/(1 - 0.5*mu**2)

    def CTbem(lambda_i):
        return(0.25*cla*sigma*((2*theta_0*(1 + (3*mu**2)/(2)))/(3) - (lambda_c + lambda_i)))
    def CTGlau(lambda_i):
        return(2*lambda_i*np.sqrt(((V*np.cos(alpha_c-a1))/(Omega*R))**2 + ((V*np.sin(alpha_c-a1))/(Omega*R) + lambda_i )**2))
    def F(lambda_i):
        return CTbem(lambda_i) - CTGlau(lambda_i)

    lambda_i = fsolve(F,[0])[0]

    a1 = ((8*mu*theta_0)/(3) - 2*mu*(lambda_c + lambda_i) - (16*q)/(gamma*Omega))/(1 - 0.5*mu**2)
    CT = CTGlau(lambda_i)

    T = CT*rho*(Omega*R)**2*np.pi*R**2
    D = CD*S*0.5*rho*V**2
    
    udot = -g*np.sin(theta_f) - (D*u)/(m*V) + (T*np.sin(theta_c-a1))/(m) - q*w
    wdot = g*np.cos(theta_f) - (D*w)/(m*V) - (T*np.cos(theta_c-a1))/(m) + q*u
    qdot = (-T*h*np.sin(theta_c-a1))/(Iy)
    theta_fdot = q
    xdot = u*np.cos(theta_f) + w*np.sin(theta_f)
    ydot = u*np.sin(theta_f) - w*np.cos(theta_f)

    u += udot*dt
    w += wdot*dt
    q += qdot*dt
    theta_f += theta_fdot*dt
    x += xdot*dt
    y += ydot*dt
    t += dt

    ulst += [u]
    wlst += [w]
    qlst += [q]
    theta_flst += [theta_f]
    tlst += [t]
    xlst += [-x]
    ylst += [y]
    ydotlst += [ydot]
    #thetaclst += [theta_c]

#plt.plot(tlst,thetaclst)
#plt.ylabel('theta_c(m/s)',rotation=0)
#plt.xlabel('t(s)')
#plt.legend
#plt.show()
plt.figure(1)
plt.plot(tlst,ulst)
plt.ylabel('u(m/s)',rotation=0)
plt.xlabel('t(s)')
plt.legend
plt.figure(2)
plt.plot(tlst,xlst)
plt.ylabel('x(m)',rotation=0)
plt.xlabel('t(s)')
plt.legend
plt.show()
plt.figure(3)
plt.plot(tlst,wlst)
plt.ylabel('w(m/s)',rotation=0)
plt.xlabel('t(s)')
plt.legend
plt.figure(4)
plt.plot(tlst,ylst)
plt.ylabel('y(m)',rotation=0)
plt.xlabel('t(s)')
plt.legend
plt.show()
plt.figure(5)
plt.plot(tlst,qlst)
plt.ylabel('q(rad/s)',rotation=0)
plt.xlabel('t(s)')
plt.legend
plt.figure(6)
plt.plot(tlst,theta_flst)
plt.ylabel(r'$\theta_f$(rad)',rotation=0)
plt.xlabel('t(s)')
plt.legend
plt.show()

