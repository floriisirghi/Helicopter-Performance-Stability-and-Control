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
    print(V)

    #---------------- Profile drag power --------------------
    mu = V/tip_speed_main #tip speed ratio

    P_p_and_P_d = (sigma*CDp)/8*rho_SL*(tip_speed_main)**3*m.pi*(R_main**2)*(1+4.65*(mu**2)) 

    #--------------- Induced power ---------------------------
    V_bar = V/v_i_hov
    #v_i_bar = 1/V_bar #works only for high speeds
    v_i_bar = np.sqrt(-((V_bar ** 2) / 2) + np.sqrt((V_bar ** 4) / 4 + 1))

    P_i = k*W*v_i_bar*np.sqrt(W/(2*rho_SL*np.pi*R_main**2)) #from equation 35a

    #----------- Total main rotor power in forward flight ---------
    P_tot_main = P_par + P_p_and_P_d + P_i



    # 2. Calculate the tail rotor power using BEM theory
    # ----------------- TAIL ROTOR-----------------------------

    # ---------------- Profile drag power --------------------
    mu_tail = V / tip_speed_tail  # tip speed ratio

    P_p_and_P_d_tail = (sigma_tail * CDp) / 8 * rho_SL * (tip_speed_tail ** 3) * m.pi * (R_tail ** 2) * (1 + 4.65 * (mu_tail ** 2))

    # --------------- Induced power ---------------------------
    V_bar = V / v_i_hov

    #v_i_bar = 1/V_bar #works only for high speeds
    v_i_bar = np.sqrt(-((V_bar ** 2) / 2) + np.sqrt((V_bar ** 4) / 4 + 1))

    T_tail = P_hov_BEM/((tip_speed_tail/R_tail)*l_tr)

    P_i_tail = 1.1*k_tail*T_tail*v_i_bar*v_i_hov  # from equation 35a

    # ----------- Total main rotor power in forward flight ---------
    P_tot_tail = P_p_and_P_d_tail + P_i_tail

    P_tot = P_tot_main + P_tot_tail

    return P_tot


V = np.linspace(0,100,101)
Power = rotor_power_forward_flight(V)/1000
plt.plot(V,Power)
plt.plot(V,0.98*1270*np.ones(len(V)),label='Max TO')
plt.plot(V,0.98*1194*np.ones(len(V)),label='Max CO')
plt.legend()
plt.title(" Power velocity curve with Max Take-Off and Max Continuous power.")
plt.show()