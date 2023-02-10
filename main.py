## -*- coding: utf-8 -*-

import rayTracing as rt_line
import radPat as rp
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import input as I
import pandas as pd
import rayTubes as rtube

h2 = I.h2
p = I.p
n2 = I.n_diec #dielectric refractive index
n1 = I.n1 #air refractive indeix 
nML = I.nML
N = I.N
L = I.L
Array = I.Array
D = I.D
m_max = I.m_max
k0 = I.k0

s1 = I.s1
s2 = I.s2
s0 = I.s0
matchingLayer1 = I.matchingLayer1
matchingLayer2 = I.matchingLayer2
surface1_arr = I.surface1
MLayer1_arr = I.MLayer1
surface2_arr = I.surface2
MLayer2_arr = I.MLayer2
aperture_plane = I.aperture_plane

class LineSegment:
    def __init__(self, normal0, normal1, n1, n2, A, B, u, t, isLast, isFirst):
        self.normal0 = normal0 #normal of the starting point of the segment
        self.normal1 = normal1 #normal of the end point of the segment        self.n1 = n1
        self.n2 = n2
        self.A = A
        self.B = B
        self.u = u
        self.t = t
        self.isLast = isLast
        self.isFirst = isFirst

class Ray:
    def __init__(self, Pk, sk, normals, ray_lengths, idxs):
        self.Pk = Pk #all the intersection points between the initial ray and the segments
        self.sk = sk #value of the last pointing vector, direction of ray coming out from the dome
        self.normals = normals #array with all the normals for each intersection point
        self.ray_lengths = ray_lengths #all distances travelled by the ray, including from  the dome to the aperture plane
        self.idxs = idxs #indexes of all intersected segments        

segments =[]
segments = np.append(segments, rt_line.discretize_function(s0, 1, 1,1, False, True))
segments = np.append(segments, rt_line.discretize_function(aperture_plane, 1, 1, 1, True, False))

if I.matchingLayers:
    segments = np.append(segments, rt_line.discretize_function(s1, 30, nML, n2, False, False))
    segments = np.append(segments, rt_line.discretize_function(matchingLayer1, 30, n1, nML, False, False))
    segments = np.append(segments, rt_line.discretize_function(matchingLayer2, 30, nML, n1, False, False))
    segments = np.append(segments, rt_line.discretize_function(s2, 30, n2, nML, False, False))
else:
    segments = np.append(segments, rt_line.discretize_function(s1, 30, n1, n2, False, False))
    segments = np.append(segments, rt_line.discretize_function(s2, 30, n2, n1, False, False))


# if 0:
#     #if we want to import an aritrary shape from a file
#     surface1 = np.loadtxt('surface1.csv', delimiter=',')
#     surface2 = np.loadtxt('surface2.csv', delimiter=',')


################################# REVERSE ##########################################3
angle_in = []
angle_position = []
angle_in, angle_position = rt_line.reverseRayTracing_segments(I.output_angle, segments)
f = interp1d(angle_position, angle_in, kind='cubic')
angles_for_direct = f(Array)
################################# END REVERSE ##########################################3



################################# DIRECT ##########################################3
angle_in = np.ones(N)*np.deg2rad(I.output_angle)
fig = plt.figure(2)
fig.set_dpi(300)
ax = fig.add_subplot(111)
csfont = {'fontname':'Times New Roman'}
plt.ylim([0,0.750])
plt.xlim([-0.8, 0.8])
plt.ylabel('z (mm)')
plt.xlabel('x (mm)')
plt.title('Direct Ray Tracing', **csfont)
plt.plot(p, surface1_arr, color='grey', linewidth = 0.5)
plt.plot(p, surface2_arr, color='grey', linewidth = 0.5)
plt.plot(p, MLayer1_arr, color = 'chocolate', linewidth = 0.1)
plt.plot(p, MLayer2_arr, color = 'chocolate', linewidth = 0.1)
ax.fill_between(p, MLayer1_arr, surface1_arr, color = 'orange')
ax.fill_between(p, MLayer2_arr, surface2_arr, color = 'orange')
ax.set_aspect(1, adjustable='box')
for i in range(0, len(segments)):
    plt.plot([segments[i].A[0], segments[i].B[0]], [segments[i].A[1], segments[i].B[1]], color = 'red', linewidth = 0.5)


Ak_ap = []
phi_a = np.zeros(N)
dck = []
rays = []
x_ap = np.zeros(N)
y_ap = np.zeros(N)

rays = rt_line.directRayTracing_segments(angles_for_direct, segments)
for i in range(0, len(rays)):
    Pk_np = rays[i].Pk
    for j in range(0, len(Pk_np)-3):
        if j % 2 == 0:
            plt.plot([Pk_np[j], Pk_np[j+2]], [Pk_np[j+1], Pk_np[j+3]], color='black', linewidth = 1)
            x_ap[i] = Pk_np[j]
            y_ap[i] = Pk_np[j+1]    

################################# END DIRECT, END GO ##########################################3



############################## RAY TUBE THEORY ##########################################3
path_length = np.zeros(N, dtype=np.complex_)
nk = np.zeros([N,2]) #normal of the aperture
sk = np.zeros([N,2]) #pointying vector
ts_coeff = np.ones(N, dtype=np.complex_)

#all these functions can be optimized
path_length = rtube.getPathLength(rays, segments)
nk, sk = rtube.getLastNormal(rays)
ts_coeff = rtube.getTransmissionCoef(rays, segments)
Ak_ap, dck = rtube.getAmplitude(rays, segments)

plt.grid()
plt.show()



############################## RADIATION PATTERN - KIRCHOFF ##########################################3
Etotal, theta = rp.getRadiationPattern(Ak_ap, path_length[1:N-1], nk[1:N-1], sk[1:N-1], dck, x_ap[1:N-1], y_ap[1:N-1], ts_coeff[1:N-1])
Etotal_dB = 20*np.log10(abs(Etotal))
print(max(Etotal_dB))

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
ax2.xaxis.label.set_fontsize(10)
ax2.yaxis.label.set_fontsize(10)
plt.grid()
plt.show()

##saving the radiation pattern results in an excel
df = pd.DataFrame(Etotal_dB, theta)
df.to_excel('RT_radpat_' + str(I.output_angle) + 'deg.xlsx', sheet_name='Sheet1')



    
    
    
    
    
    
    
    
    