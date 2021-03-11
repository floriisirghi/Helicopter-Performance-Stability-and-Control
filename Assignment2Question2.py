from Parameters import *
import math as m
import matplotlib.pyplot as plt

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

    return P_tot

print(rotor_power_forward_flight(150))

V = np.linspace(0,90,101)
Power = 1.1*rotor_power_forward_flight(V)/1000
plt.plot(V,Power)
plt.plot(V,1270*np.ones(len(V)),label='Max TO')
plt.plot(V,1194*np.ones(len(V)),label='Max CO')
plt.legend()
plt.show()