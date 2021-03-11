"""
Helicopter Performance, Stability and Control
Eurocopter AS365 N3 Dauphin
"""

import math as m
import numpy as np
#------------------- Global parameters  -----------------------
#------------------- Physical constants -----------------------
g = 9.80665         # [kgm/s^2] Gravity
rho_SL = 1.225       #[kg/m^3] density at sea level

#------------------- Helicopter parameters -----------------------
mass = 4300          #[kg] "maximum all up weight"
W = 4300*g           #[N] the actual weight of the helicopter
D_main = 11.94       #[m] main rotor diameter
R_main = D_main/2    #[m] main rotor radius
D_tail = 1.1         #[m] tail rotor diameter
P_TO = 635           #[kW] take-off power
P_cont = 597         #[kW] continuous power
Disc_L = W/(m.pi * (R_main**2)) #[N/m^2] disc loading
v_i_hov = 12.3983    #[m/s] rotor induced velocity in hover
sum_CD_S = 1         #[m^2] "flat plate area"
tip_speed_main = 218 #[m/s] this value is taken from the brightspace database, so for the older version of the helicopter
tip_speed_tail = 227 #[m/s] this value is taken from the brightspace database, so for the older version of the helicopter
N_Blades = 4         #[/] Number of blades.
c = 0.3              #[m] Average chord.                    Assumption
M = 0.7              #[/] Figure of merit.                  Assumption
CDp = 0.04           #[/] Average drag coefficient.         Assumption
Omega = tip_speed_main/(R_main)     #[/s] Rotational rate.      Could be angular, then correct with 2pi, have to investigate
sigma = N_Blades*c/(np.pi*R_main)   # Rotor solidity.       Uses assumptions as input
k = 1.15             # induced drag power factor; this parameter has a value between 1.1 and 1.2 as stated in the reader, so we can assume it to be 1.15

