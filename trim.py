import numpy as np

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

V_min = 30
V_max = 80      #Determine from assignment 2!

V_array = np.linspace(V_min,V_max,5*V_max)
theta_c = np.zeros((1,len(V_array)))
theta_0 = np.zeros((1,len(V_array)))

lambda_i = 0.5  #Initial guess
a1 = 0.5        #Initial guess
theta_0 = 0.5   #Initial guess
for i in range(len(V_array)):
    V = V_array[i]
    mu = V/(Omega*R)
    D = 0.5*rho*CDS*V**2
    W = m*g

    Ct = 0.25*cla*sigma*(2*theta_0*(1 + 1.5*mu**2)/3 - (mu*a1 + mu*(D/W) + lambda_i))
    A = np.matrix([[1+1.5*mu**2, -8*mu/3],[-mu, 2/3 + mu**2]])
    b = np.array([[-2*mu**2*(D/W) - 2*mu*lambda_i],[((4*Ct)/(sigma*cla)) + mu*(D/W) + lambda_i]])
    print(A)
    print(b)
    x = np.linalg.solve(A,b)
    print(x[0][0])
    print(x[1][0])
    theta_c[i] = x[0][0]
    theta_0[i] = x[1][0]
    a1 = theta_c[i]