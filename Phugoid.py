from Parameters import *
from MoI import Iyy_total
import numpy as np
import matplotlib.pyplot as plt
from trim import knotstomps, trimconditions, mpstoknots


tau = 0.1
#lock = (rho_SL*cl_alpha*c*R_main**4)/Iyy_total #this gives a very small lock number
lock = 6
h=1 #this is the value in the Matlab file but maybe we should change it; it's called mast in there


# ---------- Initial values --------------
t0=0
steps=3600
time=360 #sec
step=(time-t0)/steps


V_man1 = knotstomps(90)
theta_c_gen, theta_0_gen = trimconditions(knotstomps(90))
collect=[theta_0_gen] + (steps-1)*[0]
longit=[theta_c_gen] + (steps-1)*[0]

u0=knotstomps(90)
w0=0
q0=0
pitch0=0

x0=0
lambda_i0 = np.sqrt(mass*g/(np.pi*R_main**2*2*rho_SL))/tip_speed_main #nondimensional inducd velocity

t=[t0] + (steps-1)*[0]
u=[u0] + (steps-1)*[0]
w=[w0] + (steps-1)*[0]
q=[q0] + (steps-1)*[0]
pitch=[pitch0] + (steps-1)*[0]

x=[x0] + (steps-1)*[0]
lambda_i=[lambda_i0] + (steps-1)*[0]
z=[-100] + (steps-1)*[0]
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
dc = steps*[0]
V = steps*[0]
dV = steps*[0]
c_des=0
h_des = 100
pitch_des = 0
dtf = steps*[0]
altitude_h = steps*[0]

#Gains needed for the cyclic controller
K1 = 0.031
K2 = 0.12
K3 = 0.002
K4 = -0.006
K5 = 0.07
K6 = -0.0001

#Gains needed for the collective controller
K7 = 0.1
K8 = 0.05
K9 = 0.04
K10 = 0.89

V_des = knotstomps(90)
#----------------- Integration scheme -------------------
PilotOn = True
for i in range(steps):
    PilotOn = True
    phugoid = True
    if phugoid == True:
        print(t[i])
        if 299 <= t[i] <= 300:
            pitch_use = pitch[i]

        if 300 <= t[i] <=301:
            pitch[i] = pitch_use + 1*np.pi/180
            print('hello2')
    
    if PilotOn == True:

    # Law for cyclic
        V[i] = np.sqrt(u[i] ** 2 + w[i] ** 2)
        if (i + 1) < steps:
            dV[i + 1] = dV[i] + -(V[i] - V_des) * step

            pitch_des = K4 * -(V[i] - V_des) + K5 * udot[i] + K6 * dV[i]

        if (i + 1) < steps:
            dtf[i + 1] = dtf[i] + (pitch[i] - pitch_des) * step

        if i>=1:
            longit[i] = theta_c_gen + K1 * (pitch[i] - pitch_des) * 180 / np.pi + K2 * q[i] * 180 / np.pi + K3 * dtf[i]*180/np.pi

        # Law for collective
        c[i] = u[i] * np.sin(pitch[i]) - w[i] * np.cos(pitch[i])
        altitude_h[i] = -z[i]
        c_des = K9 * (h_des - altitude_h[i]) + K10 * c[i]
        if (i + 1) < steps:
            dc[i + 1] = dc[i] + (c_des - c[i]) * step
            collect[i] = theta_0_gen + K7 * (c_des - c[i]) + K8 * dc[i + 1]

    if PilotOn == False:
        print("False")
        # Law for cyclic
        V[i] = np.sqrt(u[i] ** 2 + w[i] ** 2)
        
        # Law for collective
        c[i] = u[i] * np.sin(pitch[i]) - w[i] * np.cos(pitch[i])
        altitude_h[i] = -z[i]
        c_des = K9 * (h_des - altitude_h[i]) + K10 * c[i]

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

    udot[i]= -g*np.sin(pitch[i]) - sum_CD_S/mass*0.5*rho_SL*u[i]*vv[i] + thrust[i]/mass*np.sin(helling[i]) - q[i]*w[i]

    wdot[i]= g*np.cos(pitch[i]) - sum_CD_S/mass*0.5*rho_SL*w[i]*vv[i] - thrust[i]/mass*np.cos(helling[i] + q[i]*u[i])

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


plotting = True

if plotting == True:
    plt.figure(1)
    plt.plot(t,1.94*np.ones(len(u))*u)      # This is now in knots
    plt.ylabel('u(kts)',rotation=0)
    plt.xlabel('t(s)')
    plt.legend
    #plt.figure(2)
    #plt.plot(t,x)
    #plt.ylabel('x(m)',rotation=0)
    #plt.xlabel('t(s)')
    #plt.legend
    plt.figure(3)
    plt.plot(t,w)
    plt.ylabel('w(kts)',rotation=0)     # This is now in knots
    plt.xlabel('t(s)')
    plt.legend
    plt.figure(4)
    plt.plot(t,altitude_h)
    plt.ylabel('h(m)',rotation=0) #-z
    plt.xlabel('t(s)')
    plt.legend
    #plt.figure(5)
    #plt.plot(t,q)
    #plt.ylabel('q(rad/s)',rotation=0)
    #plt.xlabel('t(s)')
    #plt.legend
    plt.figure(6)
    plt.plot(t, pitch)
    plt.ylabel('pitch (rad)',rotation=0)
    plt.xlabel('t(s)')
    plt.legend

    plt.figure(7)
    plt.plot(t, longit)
    plt.ylabel('longitudinal control (cyclic)',rotation=0)
    plt.xlabel('t(s)')
    plt.legend
    plt.show()
