import numpy as np
import matplotlib.pyplot as plt
from Parameters import *
from Assignment2Question2 import rotor_power_forward_flight

V = np.linspace(0,120,101)
def P(V):
    return (V-2)**2 + 6

def dPdV(P,V):
    dPdV = np.zeros(len(P)-1)
    for i in (range(len(P)-1)):
        dPdV[i] = (P[i+1]-P[i])/(V[i+1]-V[i])
    return dPdV

P = rotor_power_forward_flight(V)
test = dPdV(P,V)

print(test)
plt.plot(V,P)
plt.plot(V[1:],test)
plt.show()