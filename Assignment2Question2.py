from Parameters import *
import numpy as np
from Assignment1 import *
import matplotlib.pyplot as plt

def rotor_power_forward_flight(V):

    # 1. Calculate the helicopter rotor power in forward flight. For this determine the parasite drag power,
    #    induced power and total profile drag power of the rotor.

    #-----------------MAIN ROTOR-----------------------------

    #---------------- Parasite drag power -------------------
    D_par = sum_CD_S*1/2*rho_SL*(V**2)
    P_par = D_par * V

    #---------------- Profile drag power --------------------
    mu = V/tip_speed_main #tip speed ratio

    P_p_and_P_d = (sigma*CDp)/8*rho_SL*(tip_speed_main**3)*m.pi*(R_main**2)*(1+4.65*(mu**2))

    #--------------- Induced power ---------------------------
    V_bar = V/v_i_hov


    v_i_bar_low = np.sqrt(-((V_bar[0:40] ** 2) / 2) + np.sqrt((V_bar[0:40] ** 4) / 4 + 1))
    #v_i_bar_med = (np.sqrt(-((V_bar[20:40] ** 2) / 2) + np.sqrt((V_bar[20:40] ** 4) / 4 + 1)) + T/(2*rho_SL*np.pi*R_main**2*V[20:40])/v_i_hov)/2
    v_i_bar_high = T/(2*rho_SL*np.pi*R_main**2*V[40:])/v_i_hov


    v_i_bar = np.concatenate([v_i_bar_low, v_i_bar_high])

    #v_i_bar = np.sqrt(-((V_bar ** 2) / 2) + np.sqrt((V_bar ** 4) / 4 + 1))

    P_i = k*W*v_i_bar*np.sqrt(W/(2*rho_SL*np.pi*R_main**2)) #from equation 35a

    #----------- Total main rotor power in forward flight ---------
    P_tot_main = P_par + P_p_and_P_d + P_i



    # 2. Calculate the tail rotor power using BEM theory
    # ----------------- TAIL ROTOR-----------------------------

    # ---------------- Profile drag power --------------------
    mu_tail = V / tip_speed_tail  # tip speed ratio

    P_p_and_P_d_tail = (sigma_tail * CDp) / 8 * rho_SL * (tip_speed_tail ** 3) * m.pi * (R_tail ** 2) * (1 + 4.65 * (mu_tail ** 2))

    # --------------- Induced power ---------------------------
    T_tail = P_hov_BEM / ((tip_speed_tail / R_tail) * l_tr)

    v_i_hov_tail = np.sqrt(T_tail/(2*np.pi*rho_SL*R_tail**2))
    V_bar_tail = V / v_i_hov_tail

    #v_i_bar = 1/V_bar #works only for high speeds
    v_i_bar_tail = np.sqrt(-((V_bar_tail ** 2) / 2) + np.sqrt((V_bar_tail ** 4) / 4 + 1))

    v_i_tail = v_i_bar_tail*v_i_hov_tail

    P_i_tail = 1.1*k_tail*T_tail*v_i_bar_tail*v_i_hov_tail # from equation 35a

    # ----------- Total main rotor power in forward flight ---------
    P_tot_tail = P_p_and_P_d_tail + P_i_tail

    P_tot = P_tot_main + P_tot_tail

    return P_tot_main, P_p_and_P_d_tail, mu_tail, P_i_tail, T_tail, v_i_tail, P_tot, P_tot_tail, P_par, P_p_and_P_d, P_i, v_i_bar

eta_m = 0.95        #Mechanical efficiency
P_max_TO = 1270*eta_m
P_max_CO = 1194*eta_m

V = np.linspace(0,100,101)
P_tot_main, P_p_and_P_d_tail, mu_tail, P_i_tail, T_tail, v_i_tail, P_tot, P_tot_tail, P_par, P_p_and_P_d, P_i, v_i_bar = rotor_power_forward_flight(V)
Power = P_tot/1000  #why do we have this 1.1 here?
#plt.plot(V,P_tot_tail/1000, label="Tail rotor total power")
#plt.plot(V,Power, label="Total power")
#plt.plot(V, v_i_bar, label="Induced velocity")
#plt.plot(V,P_par/1000, label="Parasite drag power")
#plt.plot(V,P_p_and_P_d/1000, label="Total profile drag power")
#plt.plot(V, P_i/1000, label="Induced power")
#plt.xlabel("Velocity [m/s]")
#plt.ylabel("Power [kW]")
#plt.plot(V,1270*0.95*np.ones(len(V)),label='Max TO')
#plt.plot(V,1194*0.95*np.ones(len(V)),label='Max CO')
#plt.legend()
#plt.title(" Power velocity curve with Max Take-Off and Max Continuous power.")
#plt.show()

#print(P_par[60]+ P_i[60] + P_p_and_P_d[60])
#print(P_par[60])
#print(P_i[60])
#print(P_p_and_P_d[60])

print(P_p_and_P_d_tail[60])
print(mu_tail[60])
print(P_p_and_P_d_tail[60]+P_i_tail[60])
print(P_tot_tail[60])
print(P_tot_tail[60]*100/P_tot_main[60])

"""
eta_m = 0.95        #Mechanical efficiency
P_max_TO = 1270*eta_m
P_max_CO = 1194*eta_m

V = np.linspace(0,100,101)
Power = rotor_power_forward_flight(V)/1000
plt.plot(V,Power)
plt.xlabel("Velocity [m/s]")
plt.ylabel("Power [kW]")
plt.plot(V,P_max_TO*np.ones(len(V)),label='Max TO')
plt.plot(V,P_max_CO*np.ones(len(V)),label='Max CO')
plt.legend()
plt.title(" Power velocity curve with Max Take-Off and Max Continuous power.")
plt.show()


V = np.linspace(0,100,101)
Power = rotor_power_forward_flight(V)/1000
plt.plot(V,Power)
plt.plot(V,0.98*1270*np.ones(len(V)),label='Max TO')

#why do we have this 0.98 here?

plt.plot(V,0.98*1194*np.ones(len(V)),label='Max CO')
plt.legend()
plt.title(" Power velocity curve with Max Take-Off and Max Continuous power.")
plt.show()
"""
