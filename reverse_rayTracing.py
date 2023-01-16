# -*- coding: utf-8 -*-
"""
SMaria pUBILL
"""
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import openpyxl
from sympy import *
from scipy.optimize import fsolve
import input as I

c1 = I.c1
c2 = I.c2
k1 = I.k1
k2 = I.k2
h1 = I.h1
h2 = I.h2
p = I.p
er = I.er
mur = I.mur
e0 = I.e0
n2 = I.n_diec #dielectric refractive index
n1 = I.n1 #air refractive indeix 
nML = I.nML
wv = I.wv # wavelength in mm (defined in the paper)
k0 = I.k0 #propagation constant in free space
N = I.N
L = I.L
Array = I.Array
D = I.D
m_max = I.m_max


class Surface:
    def __init__(self, n1, n2, function):
        self.n1 = n1
        self.n2 = n2
        self.f = function

s1 = I.s1
s2 = I.s2
matchingLayer1 = I.matchingLayer1
matchingLayer2 = I.matchingLayer2

if I.nSurfaces == 4:
    surface1 = Surface(nML, n2, s1)
    MLayer1 = Surface(n1, nML, matchingLayer1)    
    surface2 = Surface(n2, nML, s2)
    MLayer2 = Surface(nML, n1, matchingLayer2)    
def s0(x):
    return 0   
surface0 = Surface(n1, n1, s0) 

   

#what you have to change
theta_o_y = I.output_angle
theta_o_x = np.deg2rad(90 -  theta_o_y) #incident angle with respect to x axis
spacing = 10
long_r3 = h2*5 #how long is the 3rd ray (from dome surface to aperture plane)



#=============================================================================
def getSurfacePoints(s,p):
    array = []
    spoints=[]
    index = 0
    s_aux = np.ones(len(p))*s if type(s) == float else s
    for i in range(0,len(s_aux)-1):
        if(i%spacing==0):
            spoints = np.append(spoints, s_aux[i])
            array = np.append(array, p[i])

        index += index
    return spoints, array
#=============================================================================


#=============================================================================
def snell(theta_out, ni, no):
    arg = no/ni * np.sin(theta_out)
    if abs(arg) <= 1:
        theta_ref = np.arcsin(no / ni * np.sin(theta_out))
    else:
        theta_ref = 0.
    return theta_ref
#=============================================================================

#=============================================================================
def getTheta_btw(m1, m2):
    return np.arctan(((m2-m1)/(1+m2*m1))) #get the angle between two slopes
#=============================================================================

#=============================================================================
def getAngleBtwVectors(v1, v2):
    return np.arctan2( v1[0]*v2[1] - v1[1]*v2[0], v1[0]*v2[0] + v1[1]*v2[1] )
#=============================================================================

#=============================================================================
def __comp_der(f,z):
    """
    Computation of the residue of f in z
    """
    #h = 1.E-12*z
    h = 1.E-6
    f1 = f(z-h)
    f2 = f(z+h)
    #return (f2-f1)/(2.*h)
    return np.gradient(np.array([f1, f(z), f2]), h)[1]
#=============================================================================

#=============================================================================
def findNormal(x,f):
    #def F(t): return f(hi, ci, ki, t)
    def F(t): return f(t)
    m_t = __comp_der(F, x)
    if abs(m_t) == 0:
        m_n = 1.E6
    else:
        m_n = -1./(m_t) #the slope of the normal line
    return m_n
#=============================================================================

#=============================================================================
def getTheta_i_max(  ):
    return  np.arcsin(n1/n2) #get the critical angle inside the dielectric with respect to normal 2
#=============================================================================

#=============================================================================
def getPhaseDisrt_i(d1,d2,d3): #get the phase distribution at the aperture plane
    return (d1+d3)*k0 + d2*k0*np.sqrt(er)
#=============================================================================


#==================================================
def distance(pointA, pointB):
    return (
        ((pointA[0] - pointB[0]) ** 2) +
        ((pointA[1] - pointB[1]) ** 2)
    ) ** 0.5  # fast sqrt
#==================================================


def findIntersection_2func(fun1,fun2,x0):
    # z = np.array([y -m3*(x-x1) - y1, y - h1 + (c1*pow(x,2))/(1+np.sqrt(1-(1+k1)*pow(c1,2)*pow(x,2))) ])
    def f(xy):    
        x, y =xy 
        return np.array([y-fun1(x), y-fun2(x)])
    return fsolve(f, [-0.2000, 0.2000], full_output=1)  



#=============================================================================
def findIntersectionv2(fun1, Pi, vi, all):
    # z = np.array([y -m3*(x-x1) - y1, y - h1 + (c1*pow(x,2))/(1+np.sqrt(1-(1+k1)*pow(c1,2)*pow(x,2))) ])
    def f(xy):    
        x, y =xy 
        return np.array([y-fun1(x), y-s1(x)])
    def g(xy):    
        x, y =xy 
        return np.array([y-fun1(x), y-s2(x)]) 
    def h(xy):    
        x, y =xy 
        return np.array([y-fun1(x), y-r3_ort(x)]) 
    def ml1(xy):    
        x, y =xy 
        return np.array([y-fun1(x), y-matchingLayer1(x)])    
    def ml2(xy):    
        x, y =xy 
        return np.array([y-fun1(x), y-matchingLayer2(x)])
    def arr(xy):
        x,y = xy
        return np.array([y-fun1(x), y-s0(x)])            

    result_s1 = fsolve(f, [-0.200, 0.200], full_output=1)
    result_s2 = fsolve(g, [-0.200, 0.200], full_output=1)
    result_s3 = fsolve(h, [-0.200, 0.200], full_output=1)
    result_ml1 = fsolve(ml1, [-0.200, 0.200], full_output=1)
    result_ml2 = fsolve(ml2, [-0.200, 0.200], full_output=1)
    result_s0 = fsolve(arr, [-0.200, 0.200], full_output=1)
  
    intersection = [[result_s0[0][0], s0(result_s0[0][0])], [result_ml1[0][0], matchingLayer1(result_ml1[0][0])], [result_s1[0][0], s1(result_s1[0][0])], [result_s2[0][0], s2(result_s2[0][0])], \
         [result_ml2[0][0], matchingLayer2(result_ml2[0][0])], [result_s3[0][0], r3_ort(result_s3[0][0])]]
    v = [[result_s0[0][0]-Pi[0], s0(result_s0[0][0])-Pi[1]], [result_ml1[0][0]-Pi[0], matchingLayer1(result_ml1[0][0])-Pi[1]], [result_s1[0][0]-Pi[0], s1(result_s1[0][0])-Pi[1]],[result_s2[0][0]-Pi[0], s2(result_s2[0][0])-Pi[1]] , \
        [result_ml2[0][0]-Pi[0], matchingLayer2(result_ml2[0][0])-Pi[1]], [result_s3[0][0]-Pi[0], r3_ort(result_s3[0][0])-Pi[1]]]     
    j = -1
    dist = 1e5
    for i in range(0, len(v)):
        #origin = np.array([Pi[0], Pi[1]])
        #plt.quiver(*origin, *v[i], color='red')
        aux = np.dot(vi, v[i])
        if dist > aux and intersection[i] != Pi and aux > 0.001: 
            dist=aux
            j = i
    
  #  sk = np.append(sk, v[j])
    if all == 1:
        return intersection
    else:
        if j == 0:
            return intersection[j], surface0
        elif j == 1:
            return intersection[j], MLayer1
        elif j == 2:
            return intersection[j], surface1
        elif j == 3:    
            return intersection[j], surface2  
        elif j == 4:
            return intersection[j], MLayer2             
        elif j == 5: 
            return intersection[j], aperture_plane       
        else: 
            return [0,0], surface0                      
#=======================================================================



#=======================================================================
def ray_reverse(vi, ri, x1, y1, iterations, nki, ski, Pki):
    iterations = iterations + 1
    [xi, yi], solution =  findIntersectionv2(ri, [x1, y1], vi, 0) 
    f = solution.f
    n_in = solution.n1
    n_out = solution.n2

       #plt.plot(xi, yi, 'bx')

    origin = np.array([x1, y1])
    m_n = findNormal(xi, f)
    v_n = -np.array([1,m_n])*np.sign(m_n) #normal vector
    v_n_norm = v_n/np.sqrt(v_n[0]**2 + v_n[1]**2) #normal unit vector
     
    origin = np.array([xi, yi])
    theta_i = getAngleBtwVectors(v_n, vi)
    theta_t = snell(theta_i, n_in, n_out)
    if [xi, yi] == [0,0] or (abs(theta_i) > getTheta_i_max() and n_out == n2): 
        return -1, [0,0], [0,0]


   # plt.quiver(*origin, *v_n, color='b')
    u = np.cos(theta_t)*v_n_norm[0] - np.sin(theta_t)*v_n_norm[1]
    v = np.sin(theta_t)*v_n_norm[0] + np.cos(theta_t)*v_n_norm[1]
    v_t = np.array([u, v])
    
    Pki= np.append(Pki,[xi, yi])

    def r_t(x):
            return (v_t[1]/v_t[0])*(x-xi)+yi
        #plt.plot(p, r_r(p), color='green', linewidth = 0.5)
        #return nk
    if iterations < I.MAX_ITERATIONS:
        nki = v_n_norm
        ski = v_t
        #all_normalsi = np.append(all_normalsi, m_n)
        #incident_anglei = np.append(incident_anglei, theta_i)

        return ray_reverse(v_t, r_t, xi, yi, iterations, nki, ski, Pki)   
    else: 
        origin = np.array([x1, y1])
    #    plt.quiver(*origin, *ski, color='b')
     #   plt.quiver(*origin, *nki, color='g')
        return Pki, ski, nki          
#=======================================================================



surface1_arr = I.surface1
MLayer1_arr = I.MLayer1
surface2_arr = I.surface2
MLayer2_arr = I.MLayer2

[surface_points, p_points] = getSurfacePoints(surface2_arr, p)
fig = plt.figure(1)
fig.set_dpi(300)
ax = fig.add_subplot(111)
plt.rcParams['font.size'] = '12'
plt.plot(p, surface1_arr, color='grey', linewidth = 0.5)
plt.plot(p, surface2_arr, color='grey', linewidth = 0.5)
ax.fill_between(p, surface1_arr, surface2_arr, color = 'lightgrey')
plt.plot(p, MLayer1_arr, color = 'chocolate', linewidth = 0.1)
plt.plot(p, MLayer2_arr, color = 'chocolate', linewidth = 0.1)
ax.fill_between(p, MLayer1_arr, surface1_arr, color = 'orange')
ax.fill_between(p, MLayer2_arr, surface2_arr, color = 'orange')
ax.set_aspect(1, adjustable='box')
plt.ylim([0,1])
plt.title('Reverse Ray Tracing')
plt.xlim([-D, D])
plt.ylabel('z (mm)' )
plt.xlabel('x (mm)')
plt.rcParams["font.family"] = "Times New Roman"

y_r_max = np.sin(theta_o_x)*h2*2
x_r_max = np.cos(theta_o_x)*D*np.sign(theta_o_y)
m3 = np.tan(theta_o_x)
if m3 > m_max: m3=m_max
m_t = -1./m3
def r3_ort(x):
    return m_t*(x - x_r_max) + y_r_max
aperture_plane = Surface(n1, n1, r3_ort)
#plt.plot(p_points, surface_points, 'bx')


row = []
Pk = [list(row) for i in range( 0, len(p_points))]
nk = np.zeros([len(p_points),2])
sk = np.zeros([len(p_points),2])
angle_in = []
phi_array = []

for i in range(0, len(p_points)): 
    #points at the aperture plane, perpendicular to ray3 
    xi_3 = p_points[i]
    yi_3 = r3_ort(xi_3)
    origin3 =    np.array([xi_3, yi_3])   
    v3 =np.array([-np.cos(theta_o_x), -np.sin(theta_o_x)])
    def r3(x):
        return (v3[1]/v3[0])*(x-xi_3)+r3_ort(xi_3)

    iterations = 0
    Pk[i] = np.append(Pk[i], [xi_3, yi_3])
    Pk[i], sk[i], nk[i] = ray_reverse(v3, r3, xi_3, yi_3, iterations, nk[i], sk[i], Pk[i])
    if isinstance(Pk[i], int): 
        continue

    x0 = Pk[i][len(Pk[i])-2]
    y0 = Pk[i][len(Pk[i])-1]
    Pk_np = np.array(Pk[i])
    if abs(x0) <= max(Array)+0.04 and y0 >= 0 :
        angle_in = np.append(angle_in, getAngleBtwVectors(sk[i], [0, -1]))
        phi_array = np.append(phi_array, x0)
        for m in range(0,I.MAX_ITERATIONS):
            plt.plot([Pk_np[m*2], Pk_np[m*2+2]], [Pk_np[m*2+1], Pk_np[m*2+3]], color='orange', linewidth = 0.5)

plt.grid()
plt.show()

#save the input angles for Direct RT to excel
df = pd.DataFrame(np.rad2deg(angle_in), phi_array)
df.to_excel('Reverse_anglesIn_'+ str(theta_o_y)+'.xlsx', sheet_name='Sheet1')

# #save the phase ditribution to excel
# df2 = pd.DataFrame(phi_a, phi_array)
# df2.to_excel('Phase_distribution_'+ str(theta_o_y) +'.xlsx', sheet_name='Sheet1')



