from Parameters import *
import math as m

def power_calcualtions(V):
# Calculate the helicopter rotor power in forward flight. For this determine the parasite drag power,
# induced power and total profile drag power of the rotor.

    #---------------- Parasite drag power -------------------
    D_par = sum_CD_S*1/2*rho_SL*(V**2)
    P_par = D_par * V

    #---------------- Profile drag power --------------------
    mu = V/omegaR_main #tip speed ratio

    P_p_and_P_d = (sigma*C_D_p)/8*rho_SL*(omegaR_main**3)*m.pi*(R_main**2)*(1+3*(mu**2)) #here we could have 4.65 instead of the 3 in front of mu

    #--------------- Induced power ---------------------------
    P_i = 