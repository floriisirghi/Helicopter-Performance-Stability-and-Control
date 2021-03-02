"""
Helicopter Performance, Stability and Control
Eurocopter AS365 N3 Dauphin
"""

import math as m

#------------------- Global parameters -----------------------
# These are constants that are defined here globally and may be imported and used when needed

g = 9.80665          #[kg/m*s^2] gravitational acceleration
mass = 4300          #[kg] "maximum all up weight"
W = 4300*g           #[N] the actual weight of the helicopter
D_main = 11.94       #[m] main rotor diameter
R_main = D_main/2    #[m] main rotor radius
D_tail = 1.1         #[m] tail rotor diameter
rho_SL = 1.225       #[kg/m^3] density at sea level
P_TO = 635           #[kW] take-off power
P_cont = 597         #[kW] continuous power
Disc_L = W/(m.pi * (R_main**2)) #[N/m^2] disc loading
v_i_hov = 12.3983    #[m/s] rotor induced velocity in hover
sum_CD_S = 1         #[m^2] "flat plate area"
tip_speed_main = 218 #[m/s] this value is taken from the brightspace database, so for the older version of the helicopter
tip_speed_tail = 227 #[m/s] this value is taken from the brightspace database, so for the older version of the helicopter