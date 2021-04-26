import numpy as np
import matplotlib.pyplot as plt
from Parameters import *

"""
Ideal Power
"""
T = W
v_i = np.sqrt(T/(2*rho_SL*np.pi*R_main**2))
P_i = T*v_i

"""
Hover Power ACT
"""
P_hov_ACT = (W/M)*np.sqrt(W/(2*np.pi*rho_SL*R_main**2))

"""
Hover Power BEM
"""
P_hov_BEM = T*v_i + (CDp/8)*rho_SL*sigma*(Omega*R_main)**3*np.pi*R_main**2

print(P_i,P_hov_ACT,P_hov_BEM)
print(P_hov_BEM/P_hov_ACT)
# 28% Difference is quite a lot, we might have to tweak the assumed parameters. Probably CDp is assumed too high.