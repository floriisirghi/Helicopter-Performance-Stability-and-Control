from Parameters import *
from MoI import Iyy_total
import numpy as np
import matplotlib.pyplot as plt

tau = 0.1
#lock = (rho_SL*cl_alpha*c*R_main**4)/Iyy_total #this gives a very small lock number
lock = 6
h=1 #this is the value in the Matlab file but maybe we should change it; it's called mast in there


# ---------- Initial values --------------
t0=0
steps=800
time=80 #sec
step=(time-t0)/steps

collect=[6*np.pi/180] + (steps-1)*[0]
longit=[0*np.pi/180] + (steps-1)*[0]

u0=0
w0=0
q0=0
pitch0=0*np.pi/180

x0=0
lambda_i0 = np.sqrt(mass*g/(np.pi*R_main**2*2*rho_SL))/tip_speed_main #nondimensional inducd velocity

t=[t0] + (steps-1)*[0]
u=[u0] + (steps-1)*[0]
w=[w0] + (steps-1)*[0]
q=[q0] + (steps-1)*[0]
pitch=[pitch0] + (steps-1)*[0]

x=[x0] + (steps-1)*[0]
lambda_i=[lambda_i0] + (steps-1)*[0]
z=[0] + (steps-1)*[0]
longitgrad=[0] + (steps-1)*[0]
q_dimless=[0] + (steps-1)*[0]
v_dimless=[0] + (steps-1)*[0]
phi=[0] + (steps-1)*[0]
alpha_c = [0] + (steps-1)*[0]
mu = [0] + (steps-1)*[0]
lambda_c = [0] + (steps-1)*[0]
a1 = [0] + (steps-1)*[0]
CT_elem = [0] + (steps-1)*[0]
alpha_d = [0] + (steps-1)*[0]
CT_glau = [0] + (steps-1)*[0]
lambda_i_dot = [0] + (steps-1)*[0]
thrust = [0] + (steps-1)*[0]
helling = [0] + (steps-1)*[0]
vv = [0] + (steps-1)*[0]

udot = [0] + (steps-1)*[0]
wdot = [0] + (steps-1)*[0]
qdot = [0] + (steps-1)*[0]
pitchdot = [0] + (steps-1)*[0]

xdot = [0] + (steps-1)*[0]
zdot = [0] + (steps-1)*[0]

c = steps*[0]
altitude_h = steps*[0]

#----------------- Integration scheme -------------------

for i in range(steps):
    if t[i]>=0.5 and t[i]<=1:
        longit[i]=1*np.pi/180
    else:
        longit[i]=0*np.pi/180

    if t[i]>=15:
        longitgrad[i]=0.2*pitch[i]*180/np.pi + 0.2*q[i]*180/np.pi #PD controller
        longit[i] = longitgrad[i]*np.pi/180

    #No law for collective

    c[i] = u[i]*np.sin(pitch[i]) - w[i]*np.cos(pitch[i])
    altitude_h[i] = -z[i]
    collect[i]=collect[0]

    # ---- Defining the differential equations ------

    # ---- Defining the nondimesnional notations -----

    q_dimless[i] = q[i]/Omega
    v_dimless[i] = np.sqrt(u[i]**2 + w[i]**2)/tip_speed_main

    if u[i]==0:
        if w[i]>0:
            phi[i]=np.pi/2
        else:
            phi[i]=-np.pi/2
    else:
        phi[i]=np.arctan(w[i]/u[i])

    if u[i]<0:
        phi[i]=phi[i]+ np.pi

    alpha_c[i]=longit[i] - phi[i]

    mu[i]=v_dimless[i]*np.cos(alpha_c[i])
    lambda_c[i]=v_dimless[i]*np.sin(alpha_c[i])

    # ------ a1 Flapping calculation --------
    a1[i] = (-16/lock*q_dimless[i] + 8/3*mu[i]*collect[i] - 2*mu[i]*(lambda_c[i] + lambda_i[i]))/(1-0.5*mu[i]**2)

    # ------ The thrust coefficient from blade elem theory ---
    CT_elem[i] = cl_alpha*sigma/4*(2/3*collect[i]*(1+1.5*mu[i]**2) - (lambda_c[i] + lambda_i[i]))

    # ------ The thrust coefficient from Glauert  ------
    alpha_d[i]=alpha_c[i]-a1[i]
    CT_glau[i]= 2*lambda_i[i]*np.sqrt((v_dimless[i]*np.cos(alpha_d[i]))**2 + (v_dimless[i]*np.sin(alpha_d[i])+lambda_i[i])**2)

    # -------------- Equations of motion ----------------
    lambda_i_dot[i] = CT_elem[i]
    thrust[i] = lambda_i_dot[i]*rho_SL*tip_speed_main**2*np.pi*R_main**2
    helling[i] = longit[i] - a1[i]
    vv[i] = v_dimless[i]*tip_speed_main

    udot[i]= -g*np.sin(pitch[i]) - sum_CD_S/mass*0.5*rho_SL*u[i]*vv[i]+thrust[i]/mass*np.sin(helling[i]) - q[i]*w[i]

    wdot[i]= g*np.cos(pitch[i]) -sum_CD_S/mass*0.5*rho_SL*w[i]*vv[i] - thrust[i]/mass*np.cos(helling[i] + q[i]*u[i])

    qdot[i]= -thrust[i]*h/Iyy_total*np.sin(helling[i])

    pitchdot[i] = q[i]

    xdot[i] = u[i]*np.cos(pitch[i])+w[i]*np.sin(pitch[i])

    zdot[i] = -c[i]

    lambda_i_dot[i]=(CT_elem[i] - CT_glau[i])/tau

    if (i+1) < steps:

        u[i+1] = u[i] + step*udot[i]
        w[i + 1] = w[i] + step * wdot[i]
        q[i + 1] = q[i] + step * qdot[i]
        pitch[i + 1] = pitch[i] + step * pitchdot[i]
        x[i+1] = x[i] + step*xdot[i]
        lambda_i[i + 1] = lambda_i[i] + step * lambda_i_dot[i]
        z[i + 1] = z[i] + step * zdot[i]
        t[i+1] = t[i] + step

plt.figure(1)
plt.plot(t,u)
plt.ylabel('u(m/s)',rotation=0)
plt.xlabel('t(s)')
plt.legend
plt.figure(2)
plt.plot(t,x)
plt.ylabel('x(m)',rotation=0)
plt.xlabel('t(s)')
plt.legend
plt.show()
plt.figure(3)
plt.plot(t,w)
plt.ylabel('w(m/s)',rotation=0)
plt.xlabel('t(s)')
plt.legend
plt.figure(4)
plt.plot(t,altitude_h)
plt.ylabel('h(m)',rotation=0) #-z
plt.xlabel('t(s)')
plt.legend
plt.show()
plt.figure(5)
plt.plot(t,q)
plt.ylabel('q(rad/s)',rotation=0)
plt.xlabel('t(s)')
plt.legend
plt.figure(6)
plt.plot(t, pitch)
plt.ylabel('pitch (rad)',rotation=0)
plt.xlabel('t(s)')
plt.legend
plt.show()

plt.figure(7)
longit_2_plot = [0]*len(longit)
for i in range(len(longit)):
    longit_2_plot[i] = longit[i]*180/np.pi
plt.plot(t, longit_2_plot)
plt.ylabel('longitudinal control (cyclic)',rotation=0)
plt.xlabel('t(s)')
plt.legend
plt.show()