"""
Helicopter Performance, Stability and Control
Eurocopter AS365 N3 Dauphin
"""

import math as m
import numpy as np
#------------------- Global parameters  -----------------------
#------------------- Physical constants -----------------------
g = 9.80665         # [kgm/s^2] Gravity
rho_SL = 1.225       #[kg/m^3] density at usual operation altitude

#------------------- Helicopter parameters -----------------------
mass = 4300          #[kg] "maximum all up weight"
W = 4300*g           #[N] the actual weight of the helicopter
D_main = 11.94       #[m] main rotor diameter
R_main = D_main/2    #[m] main rotor radius
D_tail = 1.1         #[m] tail rotor diameter
R_tail = D_tail/2    #[m] tail rotor radius
P_TO = 635           #[kW] take-off power
P_cont = 597         #[kW] continuous power
Disc_L = W/(np.pi * (R_main**2)) #[N/m^2] disc loading
v_i_hov = np.sqrt(W/(2*rho_SL*np.pi*R_main**2))    #[m/s] rotor induced velocity in hover
sum_CD_S = 1         #[m^2] "flat plate area"
tip_speed_main = 218 #[m/s] omega*R this value is taken from the brightspace database, so for the older version of the helicopter
tip_speed_tail = 227 #[m/s] this value is taken from the brightspace database, so for the older version of the helicopter
N_Blades = 4         #[/] Number of blades main rotor.
N_blades_tail = 10   #[/] Number of blades fenestron.
#c = 0.3
c = 0.385            #[m] Average chord main rotor.                    Assumption -> brightspace database
cl_alpha = 5.7 #assumed value
c_tail = 0.094       #[m] Average chord fenestron.
M = 0.7              #[/] Figure of merit.                  Assumption
Omega = tip_speed_main/(R_main)     #[/s] Rotational rate.      Could be angular, then correct with 2pi, have to investigate
sigma = N_Blades*c/(np.pi*R_main) # Rotor solidity.       Uses assumptions as input
sigma_tail = N_blades_tail*c_tail/(np.pi*R_tail) # Tail rotor solidty. Found to be 0.49 in database but for slightly older model
C_L_mean = 6.6*(W/(rho_SL*m.pi*(R_main**2)*(tip_speed_main**2)))/sigma #Medium lift coefficient
a_SL = 343           #[m/s] Speed of sound at sea level
V_cruise = 74.7      #[m/s] Recommended cruise speed
M_t = (tip_speed_main + V_cruise) /a_SL #Tip Mach number
CDp = 0.031          #[/] Mean profile drag coefficient.         Determined graphically
k = 1.15             # induced drag power factor; this parameter has a value between 1.1 and 1.2 as stated in the reader,
                     # so we can assume it to be 1.15
k_tail = 1.4         # induced drag power factor for the tail; this parameter has a value between 1.3 and 1.5
                     # as stated in the assignment description, so we can assume it to be 1.4
l_tr = 6.98          #  [m] tail rotor length. Measured on the sketch
