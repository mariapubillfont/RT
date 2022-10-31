# -*- coding: utf-8 -*-
"""
Created on Tue May 10 15:43:04 2022

@author: maria pubill
"""
from cmath import cos, log
from configparser import MAX_INTERPOLATION_DEPTH
import rayTracing as rt
import rayTracingRecursive as rtr
import radPat as rp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d
import input as I
import pandas as pd


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
Array = I.Array
MAX_ITERATIONS = I.MAX_ITERATIONS
output_angle = I.output_angle
theta_i_y = np.zeros(N)

# we read all the input angles obtained with the reverse RT and we interpolate them 
df = pd.read_excel('Reverse_anglesIn_' + str(output_angle) + '.xlsx', sheet_name='Sheet1')
df_np = np.array(df)
thy = df_np[:,1]
thy_array = df_np[:,0] 
f = interp1d(thy_array, thy, kind='cubic')
theta_i_y = -f(Array)  #the input angles in degrees with respect to y-axis

##plot the function of the input angles
# fig = plt.figure()
# fig.set_dpi(300)
# plt.plot(thy_array, thy, '.')
# plt.plot(Array, f(Array))
# plt.title('input angles from reverse RT')
# plt.grid()
# plt.show() 


#=============================================================================
def f(hi, ci, ki, p): #defining the surface shapes as conics
    return hi + (ci*pow(p,2))/(1+np.sqrt(1-(1+ki)*pow(ci,2)*pow(p,2)))
#=============================================================================


if 1:
    surface1 = f(h1, c1, k1, p)
    #surface1 = 250*np.ones(p.size) #case of flat surface
    surface1 = np.where(surface1>0, surface1, 0.)
    surface2 = f(h2, c2, k2, p)
   # surface2 = 350*np.ones(p.size) #case of flat surface
    surface2 = np.where(surface2>0, surface2, 0.)
if 0:
    #if we want to import an aritrary shape from a file
    surface1 = np.loadtxt('surface1.csv', delimiter=',')
    surface2 = np.loadtxt('surface2.csv', delimiter=',')
    
fig = plt.figure(1)
fig.set_dpi(300)
ax1 = fig.add_subplot(111)
ax1.set_aspect(1, adjustable='box')
ax1.fill_between(p, surface1, surface2, color = 'lightgrey')
plt.ylim([-500,1000])
plt.xlim([-2000,2000])
plt.ylabel('z (mm)' )
plt.xlabel('x (mm)')
plt.rcParams["font.family"] = "Times New Roman" 
ax1.xaxis.label.set_fontsize(10)
ax1.yaxis.label.set_fontsize(10)
plt.plot(p, surface1, color='grey')
plt.plot(p, surface2, color='grey')



#variables needed for the radiation pattern
nk = np.zeros([N,2]) #normal of the aperture
sk = np.zeros([N,2]) #pointying vector
Ak_ap = [] #amplitude of each ray at the lens aperture
path_length = []
dck = []
phi_a = np.zeros(N) #phase distribution oiver the array
Etotal = []
theta = []
row = []
Pk = [list(row) for i in range( 0, N)]  #all the intersection points, each row represents one ray


#plotting all the rays, directRayTracingRec(theta_i_y) from file rayTracingRecursive.py
Pk, Ak_ap, path_length, nk, sk, dck = rtr.directRayTracingRec(theta_i_y )
Pk_np = np.array(Pk)
for i in range(0,MAX_ITERATIONS):
    plt.plot([Pk_np[:, i*2], Pk_np[:, i*2+2]], [Pk_np[:, i*2+1], Pk_np[:, i*2+3]], color='black', linewidth = 0.5)

plt.grid()
plt.show()

Etotal, theta = rp.getRadiationPattern(Ak_ap, path_length[1:N-1], nk[1:N-1], sk[1:N-1], dck, Pk_np[1:N-1, 4],Pk_np[1:N-1, 5] )
Etotal_dB = 20*np.log10(abs(Etotal))

#plot the radiation pattern
fig2 = plt.figure(3)
fig2.set_dpi(400)
ax2 = fig2.add_subplot(111)
ax2.set_aspect(1.5, adjustable='box')
plt.plot(-theta*180/np.pi+90,  20*np.log10(abs(Etotal)/max(abs(Etotal))), linewidth=1, color = 'red')
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

#saving the radiation pattern results in an excel
df = pd.DataFrame(Etotal_dB, theta)
df.to_excel('RT_radpat_' + str(output_angle) + 'deg.xlsx', sheet_name='Sheet1')

   
    
    
    
    
    
    