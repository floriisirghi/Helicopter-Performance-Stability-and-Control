import numpy as np
import matplotlib.pyplot as plt
from Parameters import *
from Assignment2Question2 import rotor_power_forward_flight

def rotor_power_forward_flight(V):
# Calculate the helicopter rotor power in forward flight. For this determine the parasite drag power,
# induced power and total profile drag power of the rotor.

    #---------------- Parasite drag power -------------------
    D_par = sum_CD_S*1/2*rho_SL*(V**2)
    P_par = D_par * V

    #---------------- Profile drag power --------------------
    mu = V/tip_speed_main #tip speed ratio

    P_p_and_P_d = (sigma*CDp)/8*rho_SL*(tip_speed_main**3)*m.pi*(R_main**2)*(1+4.65*(mu**2)) 

    #--------------- Induced power ---------------------------
    V_bar = V/v_i_hov
    v_i_bar = np.sqrt(-((V_bar**2)/2)+np.sqrt((V_bar**4)/4 + 1))
    P_i = k*W*v_i_bar*np.sqrt(W/(2*rho_SL*m.pi*R_main**2)) #this is from equation 35a. Should we just use k*T*v_i as in the slides?

    #----------- Total rotor power in forward flight ---------
    P_tot = P_par + P_p_and_P_d + P_i

    return P_tot, P_par, P_p_and_P_d,P_i

V = np.linspace(0,100,101)
P_tot,P_par,P_PandP_d, P_i = rotor_power_forward_flight(V)

plt.plot(V,P_tot,label='Total')
plt.plot(V,P_par,label='Par')
plt.plot(V,P_PandP_d,label=r'P_p & P_d')
plt.plot(V,P_i,label='Induced')
plt.legend()
plt.show()

#def dPdV(V):
#    dP_pardV = 
#    dP_P_andP_ddV = 
#    dP_idV = 
#    return
