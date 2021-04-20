import numpy as np
from scipy.optimize import fsolve


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
theta_0lst = np.linspace(1,11,10)
V = 46
alpha_c = 0.2
lambda_c = 0.2
q = 0.2
mu = ((V)/(Omega*R))*np.cos(alpha_c)
lambda_i = 0.1
for theta_0 in theta_0lst:
    a1 = ((8*mu*theta_0)/(3) - 2*mu*(lambda_c + lambda_i) - (16*q)/(gamma*Omega))/(1 - 0.5*mu**2)
    print(a1)
    def CTbem(lambda_i):
        return(0.25*cla*sigma*((2*theta_0*(1 + (3*mu**2)/(2)))/(3) - (lambda_c + lambda_i)))
    def CTGlau(lambda_i):
        return(2*lambda_i*np.sqrt(((V*np.cos(alpha_c-a1))/(Omega*R))**2 + ((V*np.sin(alpha_c-a1))/(Omega*R) + lambda_i )**2))
    def F(lambda_i):
        return CTbem(lambda_i) - CTGlau(lambda_i)

    lambda_i = fsolve(F,[0])[0]
    print(lambda_i,CTbem(lambda_i),CTGlau(lambda_i))
