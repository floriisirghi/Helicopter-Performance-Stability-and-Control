import numpy as np
from scipy.optimize import fsolve

mu = 1
theta_0 = 1
lambda_c = 1
lambda_i = 1
q = 1
gamma = 1
Omega = 1

cla = 1
sigma = 1
alpha_c = 1
V = 1
R = 1

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
