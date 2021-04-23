import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

"""
Constants
"""
Vtip = 218
R = 11.94/2
Omega = Vtip/R
rho = 1.225
CDS = 1.5
m = 4300
g = 9.80665
sigma = 0.075
cla = 5.7

V_min = 0
V_max = 100      #Determine from assignment 2!
N = 500         

V_array = np.linspace(V_min,V_max,N)
theta_c = np.zeros((len(V_array)))
theta_0 = np.zeros((len(V_array)))
lambdda_i = np.zeros((len(V_array)))

Ct_init = 0

for i in range(len(V_array)):
    V = V_array[i]
    mu = V/(Omega*R)
    D = 0.5*rho*CDS*V**2
    W = m*g
    T = np.sqrt(D**2 + W**2)

    lambda_i = T/(2*m*Omega*R)
    CT = T/(rho*np.pi*R**2*(Omega*R)**2)
   
    def F(lambda_i):
        return CT - 2*lambda_i*np.sqrt((V*np.cos(D/W)/(Omega*R))**2 + (V*np.sin(D/W)/(Omega*R) + lambda_i)**2)

    lambda_i = fsolve(F,[0])[0]

    Ct = 2*lambda_i*np.sqrt((V*np.cos(D/W)/(Omega*R))**2 + (V*np.sin(D/W)/(Omega*R) + lambda_i)**2)

    A = np.matrix([[1+1.5*mu**2, -8*mu/3],[-mu, 2/3 + mu**2]])
    b = np.array([[-2*mu**2*(D/W) - 2*mu*lambda_i],[((4*Ct)/(sigma*cla)) + mu*(D/W) + lambda_i]])

    x = np.linalg.solve(A,b)

    theta_c[i] = x[0][0]
    theta_0[i] = x[1][0]
    lambdda_i[i] = lambda_i

plotting = False
if plotting == True:
    plt.plot(V_array,lambdda_i,label='lambda')
    plt.show()

    plt.plot(V_array,180/np.pi * theta_c,label=r'$\theta_c$')
    plt.plot(V_array[1:],180/np.pi * theta_0[1:],label=r'$\theta_0$')
    #plt.plot(V_array,lambdda_i,label='lambda')
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Angle [deg]')
    plt.legend()
    plt.show()


def trimconditions(V):
    mu = V/(Omega*R)
    D = 0.5*rho*CDS*V**2
    W = m*g
    T = np.sqrt(D**2 + W**2)
    lambda_i = T/(2*m*Omega*R)
    CT = T/(rho*np.pi*R**2*(Omega*R)**2)

    def F(lambda_i):
        return CT - 2*lambda_i*np.sqrt((V*np.cos(D/W)/(Omega*R))**2 + (V*np.sin(D/W)/(Omega*R) + lambda_i)**2)

    lambda_i = fsolve(F,[0])[0]

    Ct = 2*lambda_i*np.sqrt((V*np.cos(D/W)/(Omega*R))**2 + (V*np.sin(D/W)/(Omega*R) + lambda_i)**2)

    A = np.matrix([[1+1.5*mu**2, -8*mu/3],[-mu, 2/3 + mu**2]])
    b = np.array([[-2*mu**2*(D/W) - 2*mu*lambda_i],[((4*Ct)/(sigma*cla)) + mu*(D/W) + lambda_i]])

    x = np.linalg.solve(A,b)
    
    theta_c = x[0][0]
    theta_0 = x[1][0]
    
    return theta_c, theta_0

def knotstomps(knots):
    mps = 0.514*knots
    return mps

