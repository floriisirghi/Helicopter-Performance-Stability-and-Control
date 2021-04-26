# Developed by Florina Sirghi as part of the AE4314 Helicopter Performance, Stability and Control course
# April 2021

# The helicopter will be split up into the following parts:
# - main rotor, which will be approximated as a disk 2.7%
# - tail rotor, which will also be approximated as a disk -> 1.5% weight
# - fuselage (+crew and luggage), which will be approximated as an parallelepiped + cone -> 42% parallelepiped  + 8% cone
# - fuel tank, which will be approximated as a cuboid -> density of kerosene
# is 0.783 kg/l -> we carry 1135 l -> 888 kg -> appx 20% max weight+ structure -> 25%
# - engine + nacelle -> 7 %
# - main rotor hardware -> 3.8 %
# - horizontal stabilizer, which will be approximated as a rectangular plate-> 2%
# - vertical fin, which will be approximated as a parallelepiped -> 4%
# - landing gear, also approximated as a cuboid -> 4% weight

# The position of the centre of gravity considered is at 4 m from the nose and 1.2 m from the ground



from Parameters import *
import numpy as np

y_cg = 1.3
x_cg = 4
max_weight = 4300
# ------------- Main Rotor -------------

J_main = 2090 # Polar moment of inertia from Prouty (brightspace databse)
m_main = 2*J_main/(R_main**2)
y_cg_main_rotor = 3.31 #measured on the technical drawing
x_cg_main_rotor = 4
dist_y_main_2_cg = y_cg_main_rotor - y_cg
dist_x_main_2_cg = x_cg_main_rotor - x_cg

Iy_disk_main = 1/4*m_main*R_main**2

Iy_disk_main_complete = Iy_disk_main + m_main*(dist_x_main_2_cg**2 + dist_y_main_2_cg**2)

# ------------- Tail Rotor and surrounding plate -------------

m_tail = 1.5/100*4300
y_cg_tail_rotor = 1.58 #measured on the technical drawing
x_cg_tail_rotor = 10.5
dist_y_tail_2_cg = y_cg_tail_rotor - y_cg
dist_x_tail_2_cg = x_cg_tail_rotor - x_cg

Iy_disk_tail = 1/2*m_tail*(R_tail+0.3)**2

Iy_disk_tail_complete = Iy_disk_tail + m_tail*(dist_x_tail_2_cg**2 + dist_y_tail_2_cg**2)

# -------------- Fuselage, crew & luggage parallelipiped -----------

m_fuselage_p  = 42/100*max_weight
y_cg_fuselage_p  = 1.2
x_cg_fuselage_p  = 3.08
dist_y_fuselage_p_2_cg  = y_cg_fuselage_p - y_cg
dist_x_fuselage_p_2_cg  = x_cg_fuselage_p - x_cg

Iy_fuselage_p = 1/12*m_fuselage_p*(1.7**2 + 6.4**2)

Iy_fuselage_p_complete =  Iy_fuselage_p + m_fuselage_p*(dist_x_fuselage_p_2_cg**2 + dist_y_fuselage_p_2_cg**2)

# -------------- Fuselage cone -----------

m_fuselage_c  = 8/100*max_weight
y_cg_fuselage_c  = 1.25
x_cg_fuselage_c  = 7.3
dist_y_fuselage_c_2_cg  = y_cg_fuselage_c - y_cg
dist_x_fuselage_c_2_cg  = x_cg_fuselage_c - x_cg

Iy_fuselage_c = m_fuselage_c*(3/20*0.58**2 + 3/80*2.77**2)

Iy_fuselage_c_complete =  Iy_fuselage_c + m_fuselage_c*(dist_x_fuselage_c_2_cg**2 + dist_y_fuselage_c_2_cg**2)

# --------------- Fuel tank ----------------

m_fuel_tank = 25/100*max_weight
y_cg_fuel_tank = 0.55
x_cg_fuel_tank = 3.8
dist_y_fuel_tank_2_cg = y_cg_fuel_tank - y_cg
dist_x_fuel_tank_2_cg = x_cg_fuel_tank - x_cg

Iy_fuel_tank = 1/12*m_fuselage_p*(0.45**2 + 4**2)
Iy_fuel_tank_complete = Iy_fuel_tank + m_fuel_tank*(dist_x_fuel_tank_2_cg**2 + dist_y_fuel_tank_2_cg**2)

# --------------- Engine + nacelle -----------

m_eng_nac = 7/100*max_weight
y_cg_eng_nac = 2.38
x_cg_eng_nac = 4.8
dist_y_eng_nac_2_cg = y_cg_eng_nac - y_cg
dist_x_eng_nac_2_cg = x_cg_eng_nac - x_cg

Iy_eng_nac = 1/12*m_eng_nac*(4.2**2 + 0.65**2)
Iy_eng_nac_complete = Iy_eng_nac + m_eng_nac*(dist_x_eng_nac_2_cg**2 + dist_y_eng_nac_2_cg**2)


# ---------------- Main rotor hardware ---------

m_main_hard = 3.8/100*max_weight
y_cg_main_hard = 2.8
x_cg_main_hard = 4.5
dist_y_main_hard_2_cg = y_cg_main_hard - y_cg
dist_x_main_hard_2_cg = x_cg_main_hard - x_cg

Iy_main_hard = 1/12*m_main_hard*(2.1**2 + 0.4**2)
Iy_main_hard_complete = Iy_main_hard + m_main_hard*(dist_x_main_hard_2_cg**2 + dist_y_main_hard_2_cg**2)


# ---------------- Horizontal stabiliser ----------

m_hs = 2/100*max_weight
y_cg_hs = 1.43
x_cg_hs = 9.06
dist_y_hs_2_cg = y_cg_hs - y_cg
dist_x_hs_2_cg = x_cg_hs - x_cg

Iy_hs = 1/12*m_hs*0.55**2   #width of stabiliser = 0.55 m
Iy_hs_complete = Iy_hs + m_hs*(dist_x_hs_2_cg**2 + dist_y_hs_2_cg**2)

# ----------------- Vertical fin -------------------

m_vf = 4/100*max_weight
y_cg_vf = 3.1
x_cg_vf = 10.7
dist_y_vf_2_cg = y_cg_vf - y_cg
dist_x_vf_2_cg = x_cg_vf - x_cg

Iy_vf = 1/12*m_vf*(1.08**2 + 1.48**2)
Iy_vf_complete = Iy_vf + m_vf*(dist_x_vf_2_cg**2 + dist_y_vf_2_cg**2)

# ------------------ Landing gear ---------------------

m_lg = 4/100*max_weight
y_cg_lg = 0.2
x_cg_lg = 2.9
dist_y_lg_2_cg = y_cg_lg - y_cg
dist_x_lg_2_cg = x_cg_lg - x_cg

Iy_lg = 1/12*m_lg*(3.64**2 + 0.4**2)
Iy_lg_complete = Iy_lg + m_lg*(dist_x_lg_2_cg**2 + dist_y_lg_2_cg**2)


# ------------------ Total moment of inertia of the helicopter ---------------

Iyy_total = Iy_disk_main_complete + Iy_disk_tail_complete + Iy_fuselage_p_complete + Iy_fuselage_c_complete + Iy_fuel_tank_complete + Iy_eng_nac_complete + Iy_main_hard_complete + Iy_hs_complete + Iy_vf_complete + Iy_lg_complete

verify_total_mass = max_weight - (m_main + m_tail + m_fuselage_p + m_fuselage_c + m_fuel_tank + m_eng_nac + m_main_hard + m_hs + m_vf + m_lg)

#print(Iyy_total)
#print(verify_total_mass)