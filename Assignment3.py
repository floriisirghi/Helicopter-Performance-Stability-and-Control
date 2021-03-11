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
    dP_P_andP_ddV = 2*(sigma*CDp/8)*rho_SL*tip_speed_main**3*np.pi*R_main**2*4.65*V/(tip_speed_main**2)
    c = k*W*np.sqrt(W/(2*rho_SL*np.pi*R_main**2))
    a = 1/(4*v_i_hov**4)
    b = 1/(2*v_i_hov**2)
    dP_idV = c*((2*a*V**3)/(np.sqrt(a*V**4 + 1)) - 2*b*V)/(2*np.sqrt(np.sqrt(a*V**4 + 1) - b*V**2))
    dP_totdV = dP_pardV + dP_P_andP_ddV + dP_idV
    return dP_totdV/1000, dP_pardV/1000,dP_P_andP_ddV/1000,dP_idV/1000

def dPdVtest(P,V):
    dPdV = np.ones(len(V)-1)
    for i in range(len(dPdV)):
        dPdV[i] = (P[i+1]- P[i])/(V[i+1]-V[i])
    return dPdV

V = np.linspace(0,100,1001)
P_tot,P_par,P_PandP_d, P_i = rotor_power_forward_flight(V)

plt.plot(V,P_tot,label='Total')
plt.plot(V,P_par,label='Par')
plt.plot(V,P_PandP_d,label=r'P_p & P_d')
plt.plot(V,P_i,label='Induced')
plt.legend()
plt.title("Power velocity curve with individual components.")
plt.show()

analytical = dPdV(V)[0]
discrete = dPdVtest(P_tot,V)

plt.plot(V,analytical,label='analytical')
plt.plot(V[1:],discrete,label='discrete')
plt.legend()
plt.title("Verification of the analytical derivative function.")
plt.show()

def Vrangefunc(V):
    return (rotor_power_forward_flight(V)[0])/V - dPdV(V)[0]

V_range = fsolve(Vrangefunc,[50])

print("Velocity for optimum range = ",V_range[0]," m/s.")
slope = dPdV(V_range)[0]

def V_endfunc(V):
    return dPdV(V)[0]

V_end = fsolve(V_endfunc,[40])
print("Velocity for optimum endurance = ",V_end[0]," m/s.")

plt.plot(V,P_tot)
plt.plot(V,rotor_power_forward_flight(V_end)[0]*np.ones(len(V)))
plt.plot(V,slope*V)
plt.xlabel("V [m/s]")
plt.ylabel("P [kW]",rotation=0)
plt.title("Total power plot with points for max endurance and range")
plt.show()
