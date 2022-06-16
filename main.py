# -*- coding: utf-8 -*-
"""
Created on Tue May 10 15:43:04 2022

@author: maria
"""
import rayTracing as rt
import radPat as rp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import input as I



# parameters to define the conic shapes of the dome (all parameters defined in the paper)
c1 = I.c1
c2 = I.c2
k1 = I.k1
k2 = I.k2
h1 = I.h1
h2 = I.h2


# parameters to define the Array
N = I.N
L = I.L
p = I.p
er = I.er
n2 = I.n2 #dielectric refractive index
n1 = I.n1 #air refractive indeix 
wv = I.wv # wavelength in mm (defined in the paper)
k0 = I.k0
phi_a = np.zeros(N)
Etotal = []
theta = []

theta_i_y = np.zeros(N)
for i in range(0, N): theta_i_y[i] = -60

angle_out = []
m_max = 10000


#=============================================================================
def f(hi, ci, ki, p): #defining the surface shapes as conics
    return hi + (ci*pow(p,2))/(1+np.sqrt(1-(1+ki)*pow(ci,2)*pow(p,2)))
#=============================================================================


if 1:
    surface1 = f(h1, c1, k1, p)
    surface1 = np.where(surface1>0, surface1, 0.)
    # np.savetxt('surface1.csv', surface1, delimiter=',')    
    surface2 = f(h2, c2, k2, p)
    surface2 = np.where(surface2>0, surface2, 0.)
    # np.savetxt('surface2.csv', surface2, delimiter=',')       
if 0:
    surface1 = np.loadtxt('surface1.csv', delimiter=',')
    surface2 = np.loadtxt('surface2.csv', delimiter=',')
    
fig = plt.figure(1)
fig.set_dpi(300)
ax1 = fig.add_subplot(111)
ax1.set_aspect(1, adjustable='box')
ax1.fill_between(p, surface1, surface2, color = 'lightgrey')
plt.ylim([-500,h2*3])
plt.ylabel('z (mm)' )
plt.xlabel('x (mm)')
plt.rcParams["font.family"] = "Times New Roman" 
# mpl.rcParams.update({'font.size': 8})
# mpl.rcParams.update({'font.family': "Times New Roman" })
ax1.xaxis.label.set_fontsize(10)
ax1.yaxis.label.set_fontsize(10)
plt.plot(p, surface1, color='grey')
plt.plot(p, surface2, color='grey')


plt.grid()

#variables needed for the radiation pattern
nk = np.zeros([N,2]) #normal of the aperture
sk = np.zeros([N,2]) #pointying vector
Ak_ap = []
Pk = np.zeros([N,2])
Pk_intersection1 = np.zeros([N,2])
Pk_ap = np.zeros([N,2])
Pk_final = np.zeros([N,2])
path_length = []
dck = []
theta_k = []
angle_out = []
dR= []


Pk, Pk_intersection1, Pk_ap, Pk_final, sk, nk, path_length, Ak_ap, dck, theta_k, angle_out = rt.directRayTracing(surface1, surface2, theta_i_y)

plt.plot([Pk[:,0], Pk_intersection1[:,0] ], [Pk[:,1], Pk_intersection1[:,1]], color='black', linewidth = 0.5)
plt.plot([Pk_ap[:,0], Pk_intersection1[:,0] ], [Pk_ap[:,1], Pk_intersection1[:,1]], color='black', linewidth = 0.5)
plt.plot([Pk_ap[:,0], Pk_final[:,0] ], [Pk_ap[:,1], Pk_final[:,1]], color='black', linewidth = 0.5)

plt.show()
Array = np.linspace (-L/2, L/2, N)


Etotal, theta, dR = rp.getRadiationPattern(Ak_ap, path_length[1:N-1], nk[1:N-1], sk[1:N-1], dck, Pk_ap[1:N-1])
Etotal_dB = 20*np.log10(abs(Etotal)/max(abs(Etotal)))


#plot the radiation pattern
fig2 = plt.figure(3)
fig2.set_dpi(400)
ax2 = fig2.add_subplot(111)
ax2.set_aspect(1.5, adjustable='box')

plt.plot(theta*180/np.pi-90, Etotal_dB, linewidth=1, color = 'red')
plt.ylabel('Normalized Pattern, dB')
plt.xlim([-70, 70])
plt.ylim([-35, 0])
plt.xlabel('$\u03B8 $, degrees')
plt.xticks(range(-90, 91, 10))
plt.yticks(range(-35, 10, 5))
plt.rcParams["font.family"] = "Times New Roman" 
ax1.xaxis.label.set_fontsize(10)
ax1.yaxis.label.set_fontsize(10)
mpl.rcParams.update({'font.size': 10})


plt.grid()
plt.show()
    
    
    
    
    
    
    
    
    