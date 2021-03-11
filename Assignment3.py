import numpy as np
import matplotlib.pyplot as plt
from Parameters import *

V = np.linspace(-1,5,101)
def P(V):
    return (V-2)**2 + 6

def dPdV(P,V):
    dPdV = np.zeros(len(P)-1)
    for i in (range(len(P)-1)):
        dPdV[i] = (P[i+1]-P[i])/(V[i+1]-V[i])
    return dPdV
Plst = P(V)

test = dPdV(Plst,V)

print(test)
plt.plot(V,P(V))
plt.plot(V[1:],test)
plt.show()