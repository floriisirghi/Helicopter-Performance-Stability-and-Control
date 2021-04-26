from Assignment1 import P_hov_BEM
import numpy as np
import matplotlib.pyplot as plt
from Parameters import *
from Assignment2Question2 import rotor_power_forward_flight
from scipy.optimize import fsolve
plt.rcParams.update({'font.size': 16})

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

    return P_tot/1000, P_par/1000, P_p_and_P_d/1000,P_i/1000

def dPdV(V):
    dP_pardV = (3*sum_CD_S*rho_SL*V**2)/2
    dP_P_andP_ddV = 2*(sigma*CDp/8)*rho_SL*(tip_speed_main)**3*np.pi*R_main**2*4.65*V/(tip_speed_main**2)
    c = k*W*np.sqrt(W/(2*rho_SL*np.pi*R_main**2))
    a = 1/(4*v_i_hov**4)
    b = 1/(2*v_i_hov**2)
    dP_idV = c*((2*a*V**3)/(np.sqrt(a*V**4 + 1)) - 2*b*V)/(2*np.sqrt(np.sqrt(a*V**4 + 1) - b*V**2))
    dP_totdV = dP_pardV + dP_P_andP_ddV + dP_idV
    return dP_totdV/1000, dP_pardV/1000,dP_P_andP_ddV/1000,dP_idV/1000

def tailpower(V):
    # ----------------- TAIL ROTOR-----------------------------

    # ---------------- Profile drag power --------------------
    mu_tail = V / tip_speed_tail  # tip speed ratio

    P_p_and_P_d_tail = (sigma_tail * CDp) / 8 * rho_SL * (tip_speed_tail ** 3) * m.pi * (R_tail ** 2) * (1 + 4.65 * (mu_tail ** 2))

    # --------------- Induced power ---------------------------
    V_bar = V / v_i_hov

    #v_i_bar = 1/V_bar #works only for high speeds
    v_i_bar = np.sqrt(-((V_bar ** 2) / 2) + np.sqrt((V_bar ** 4) / 4 + 1))

    T_tail = (1000*rotor_power_forward_flight(V)[0])/((tip_speed_tail/R_tail)*l_tr)

    P_i_tail = 1.1*k_tail*T_tail*v_i_bar*v_i_hov  # from equation 35a

    # ----------- Total main rotor power in forward flight ---------
    P_tot_tail = P_p_and_P_d_tail + P_i_tail

    return P_tot_tail/1000

def tailpower2(V):
    # ----------------- TAIL ROTOR-----------------------------

    # ---------------- Profile drag power --------------------
    mu_tail = V / tip_speed_tail  # tip speed ratio

    P_p_and_P_d_tail = (sigma_tail * CDp) / 8 * rho_SL * (tip_speed_tail ** 3) * m.pi * (R_tail ** 2) * (1 + 4.65 * (mu_tail ** 2))

    # --------------- Induced power ---------------------------
    V_bar = V / v_i_hov

    #v_i_bar = 1/V_bar #works only for high speeds
    v_i_bar = np.sqrt(-((V_bar ** 2) / 2) + np.sqrt((V_bar ** 4) / 4 + 1))

    T_tail = (P_hov_BEM)/((tip_speed_tail/R_tail)*l_tr)

    P_i_tail = 1.1*k_tail*T_tail*v_i_bar*v_i_hov  # from equation 35a

    # ----------- Total main rotor power in forward flight ---------
    P_tot_tail = P_p_and_P_d_tail + P_i_tail

    return P_tot_tail/1000

def dPdVtail(V):
    dPppd_taildV = (sigma_tail*CDp/8)*rho_SL*(tip_speed_tail)**3*np.pi*R_tail**2*4.65*2*V/tip_speed_tail**2
    c = 1.1*k_tail*v_i_hov*(1000*rotor_power_forward_flight(V)[0])/(l_tr*(tip_speed_tail/R_tail))
    a = 1/(4*v_i_hov**4)
    b = 1/(2*v_i_hov**2)
    dPitaildV = c*(((2*a*V**3)/(np.sqrt(a*V**4 + 1))) - 2*b*V)/(2*np.sqrt(np.sqrt(a*V**4 + 1) - b*V**2))
    dPtaildV = dPppd_taildV + dPitaildV
    return dPtaildV/1000

def dPdVtail2(V):
    dPppd_taildV = (sigma_tail*CDp/8)*rho_SL*(tip_speed_tail)**3*np.pi*R_tail**2*4.65*2*V/tip_speed_tail**2
    c = 1.1*k_tail*v_i_hov*(P_hov_BEM)/(l_tr*(tip_speed_tail/R_tail))
    a = 1/(4*v_i_hov**4)
    b = 1/(2*v_i_hov**2)
    dPitaildV = c*(((2*a*V**3)/(np.sqrt(a*V**4 + 1))) - 2*b*V)/(2*np.sqrt(np.sqrt(a*V**4 + 1) - b*V**2))
    dPtaildV = dPppd_taildV + dPitaildV
    return dPtaildV/1000

def dPdVtest(P,V):
    dPdV = np.ones(len(V)-1)
    for i in range(len(dPdV)):
        dPdV[i] = (P[i+1]- P[i])/(V[i+1]-V[i])
    return dPdV

V = np.linspace(0,100,1001)
P_tot,P_par,P_PandP_d, P_i = rotor_power_forward_flight(V)



analytical = dPdV(V)[0]
discrete = dPdVtest(P_tot,V)


tail_analytical = dPdVtail(V)
tail_discrete = dPdVtest(tailpower(V),V)



def Vrangefunc(V):
    Totalpower = rotor_power_forward_flight(V)[0] + tailpower(V)
    dTotalpowerdV = dPdV(V)[0] + dPdVtail(V)
    return (Totalpower/V - dTotalpowerdV)

V_range = fsolve(Vrangefunc,[50])

print("Velocity for optimum range = ",V_range[0]," m/s.")
slope = dPdV(V_range)[0] + dPdVtail(V)

def V_endfunc(V):
    return dPdV(V)[0] + dPdVtail(V)

V_end = fsolve(V_endfunc,[40])
print("Velocity for optimum endurance = ",V_end[0]," m/s.")

Totalpower = rotor_power_forward_flight(V)[0] + tailpower(V)

plt.plot(V,Totalpower)
plt.plot(V[200:500],np.ones(len(V))[200:500]*(rotor_power_forward_flight(V_end)[0] + tailpower(V_end)))
plt.plot(V,slope*V)
plt.scatter(35.7,774.2, s=80, facecolors='none', edgecolors='r',label='Maximum Endurance')
plt.scatter(72,1055, s=80, facecolors='none', edgecolors='m',label='Maximum Range')
plt.xlabel("V [m/s]")
plt.ylabel("P [kW]",rotation=0)
plt.grid()
plt.legend()
plt.title("Total power plot with points for max endurance and range.")
plt.show()
