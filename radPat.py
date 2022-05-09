# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:29:39 2022

@author: mapf
"""
import numpy as np
import matplotlib.pyplot as plt


theta = 0
Ak_ap = np.array([0.940926, 0.961937, 0.979777, 0.989435, 0.989432, 0.979797, 0.961961, 0.940959])
dCk = np.array([86.9274,83.0587,79.9477,78.323,78.3225,79.9463,83.0564,86.9242])
rayL = np.array([314.457, 326.682, 335.695, 341.661, 344.63, 344.631, 341.672, 335.711, 326.709 ,314.496])
nk = np.transpose(np.array([[-0.174346, -0.136268, -0.0976163, -0.0585473, -0.0194465, 0.0197467, 0.058122, 0.0973439, 0.135965, 0.173961], [0.984684, 0.990672, 0.995224, 0.998285, 0.999811, 0.999805, 0.998309, 0.995251, 0.990714, 0.984753]]))
sk = np.array([[-0.268714,	0.96322], [-0.22222,	0.974997], [-0.167206,	0.985922], [-0.104074,	0.99457], [-0.0356668,	0.999364], [0.0351997,	0.99938], [0.104173,	0.994559], [0.16725,	0.985914], [0.222298,	0.974979],[0.268814,	0.963192]])



# parameters to define the conic shapes of the dome (all parameters defined in the paper)
c1 = -0.0021
c2 = -0.0005
k1 = -1.2
k2 = -3.9
h1 = 325
h2 = 345

D = 1500.
p = np.linspace(-D, D, 10000) 
er = 2.5
n2 = np.sqrt(er) #dielectric refractive index
n1 = 1. #air refractive indeix 
wv = 23. # wavelength in mm (defined in the paper)
k0 = 2*np.pi/wv #propagation constant in free space
phi_a = np.zeros(N)
theta_i_x_arr = np.deg2rad(90-theta_i_y)
angle_out = []
m_max = 10000
fig = plt.figure(1)
fig.set_dpi(300)
ax = fig.add_subplot(111)
plt.plot(p, surface1, color='grey')
plt.plot(p, surface2, color='grey')
ax.set_aspect(1, adjustable='box')
ax.fill_between(p, surface1, surface2, color = 'lightgrey')
plt.ylim([0,h2*3])
plt.ylabel('z (mm)' )
plt.xlabel('x (mm)')
plt.rcParams["font.family"] = "Times New Roman"
#plt.plot(Array, np.zeros(N), 'o', color='black')


if 1:
    surface1 = f(h1, c1, k1, p)
    surface1 = np.where(surface1>0, surface1, 0.)
    # np.savetxt('surface1.csv', surface1, delimiter=',')
    
        
    surface2 = f(h2, c2, k2, p)
    surface2 = np.where(surface2>0, surface2, 0.)
    # np.savetxt('surface2.csv', surface2, delimiter=',')
    



#=============================================================================
def getUnitVector(x1,y1, x2, y2):
    vector=[x2-x1, y2-y1]
    norma = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    u = vector/norma
    return(u)
#=============================================================================


def getRadiationPattern(Ak_ap, rayL, nk, sk, dCk, Pk_ap):
    rk =[0, 1000]
    plt.plot(0, 1000, 'x')
    
    for i in range(0,len(Ak_ap)): 
        rk_vector = getUnitVector(Pk_ap[i+2][0], Pk_ap[i+2][1], rk[0], rk[1])
        Ek = 
        # plt.quiver(Pk_ap[i+2][0], Pk_ap[i+2][1], rk_vector[0], rk_vector[1], color=['blue'], scale=15)
        E = Ak_ap[i]
    
    
    return