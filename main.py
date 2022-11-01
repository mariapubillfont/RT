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
n1 = I.n1 #air refractive indeix 
wv = I.wv # wavelength in mm (defined in the paper)
k0 = I.k0
Array = I.Array
MAX_ITERATIONS = I.MAX_ITERATIONS
phi_a = np.zeros(N)
Etotal = []
theta = []
output_angle = I.output_angle

theta_i_y = np.ones(N)*output_angle

#x = np.linspace(-L/2, L/2, 12)

# df = pd.read_excel('Reverse_anglesIn_' + str(output_angle) + '.xlsx', sheet_name='Sheet1')
# df_np = np.array(df)
# thy = df_np[:,1]
# thy_array = df_np[:,0]
    
# # row = []
# # Pk = [ list(row) for i in range( 0, N)]    
  
# #plot the function of the phases
# #x = np.linspace(-L/2, L/2, len(thy))   
# f = interp1d(thy_array, thy, kind='cubic')
# #xnew = (np.linspace(-L/2, L/2, num=1001, endpoint=True))
# # def f(x):
# #     return output_angle +x*0

# theta_i_y = f(Array)    


# fig = plt.figure()
# fig.set_dpi(300)
# plt.plot(thy_array, thy, '.')
# plt.plot(Array, f(Array))
# plt.title('input angles from reverse RT')
# plt.grid()
# plt.show() 
#theta_i_y = thy
#for i in range(0, N): theta_i_y[i] = -20

angle_out = []
m_max = 10000


#=============================================================================
def f(hi, ci, ki, p): #defining the surface shapes as conics
    return hi + (ci*pow(p,2))/(1+np.sqrt(1-(1+ki)*pow(ci,2)*pow(p,2)))
#=============================================================================


if 1:
    #surface1 = f(h1, c1, k1, p)
    surface1 = h1*np.ones(p.size)
    surface1 = np.where(surface1>0, surface1, 0.)
    # np.savetxt('surface1.csv', surface1, delimiter=',')    
    #surface2 = f(h2, c2, k2, p)
    
    surface2 = h2*np.ones(p.size)
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
Ak_ap = []
Pk = np.zeros([N,2])

Pk_intersection1 = np.zeros([N,2])
Pk_ap = np.zeros([N,2])
Pk_final = np.zeros([N,2])
path_length = []
dck = []
theta_k = []
angle_out = []
ts_coeff = np.ones(N, dtype=np.complex_)
tp_coeff = np.ones(N)
#dR= []


Pk, Ak_ap, path_length, nk, sk, dck, ts_coeff, tp_coeff = rtr.directRayTracingRec(theta_i_y)
Pk_np = np.array(Pk)
for i in range(0,MAX_ITERATIONS):
    plt.plot([Pk_np[:, i*2], Pk_np[:, i*2+2]], [Pk_np[:, i*2+1], Pk_np[:, i*2+3]], color='black', linewidth = 0.5)

plt.grid()
plt.show()



# plt.figure(2)
# plt.plot(Array, angle_out)
# plt.xlabel('Array [mm]')
# plt.ylabel('Angle out [deg]')
# plt.grid()
# plt.show()
# df = pd.DataFrame(angle_out)
# df.to_excel('theta_out.xlsx', sheet_name='Sheet1')


Etotal, theta = rp.getRadiationPattern(Ak_ap, path_length[1:N-1], nk[1:N-1], sk[1:N-1], dck, Pk_np[1:N-1, 4], Pk_np[1:N-1, 5], ts_coeff[1:N-1])
Etotal_dB = 20*np.log10(abs(Etotal))
# # Etotal_dB = 20*np.log10(abs(Etotal)/max(abs(Etotal))) + 10*np.log10(Leff/Leff_broadside)

# Pk, Pk_intersection1, Pk_ap, Pk_final, sk, nk, path_length, Ak_ap, dck, theta_k, angle_out, phi_a, Leff = rt.directRayTracing(surface1, surface2, theta_i_y, thy_array )

# plt.plot([Pk[:,0], Pk_intersection1[:,0] ], [Pk[:,1], Pk_intersection1[:,1]], color='black', linewidth = 0.5)
# plt.plot([Pk_ap[:,0], Pk_intersection1[:,0] ], [Pk_ap[:,1], Pk_intersection1[:,1]], color='black', linewidth = 0.5)
# plt.plot([Pk_ap[:,0], Pk_final[:,0] ], [Pk_ap[:,1], Pk_final[:,1]], color='black', linewidth = 0.5)


Array = np.linspace (-L/2, L/2, N)


# Etotal, theta = rp.getRadiationPattern(Ak_ap, path_length[1:N-1], nk[1:N-1], sk[1:N-1], dck, Pk_ap[1:N-1])
# Etotal_dB = 20*np.log10(abs(Etotal))
# Etotal_dB = 20*np.log10(abs(Etotal)/max(abs(Etotal))) + 10*np.log10(Leff/Leff_broadside)


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

df = pd.DataFrame(Etotal_dB, theta)
df.to_excel('RT_radpat_' + str(output_angle) + 'deg.xlsx', sheet_name='Sheet1')


plt.grid()
plt.show()



# # fig3 = plt.figure(4)
# # fig3.set_dpi(400)
# # plt.plot(Array, phi_a)
# # plt.yticks([-80, -40, 0, 40, 80], ['-80', '-40', '0', '40', '80'])
# # plt.xticks([-L/2, -L/4, 0, L/4, L/2], ['-L/2', '-L/4', '0', 'L/4', 'L/2'])
# # plt.ylim([-90,90])
# # plt.ylabel('$\phi_a$ (rad)' )
# # plt.xlabel('x (mm)')
# # plt.rcParams["font.family"] = "Times New Roman"    
# # plt.grid()    
# # plt.title('Phase distribution over the array for '+ str(output_angle) + '°')
    
    
    
    
    
    
    
    
    