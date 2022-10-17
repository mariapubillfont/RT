# -*- coding: utf-8 -*-
"""
maria pubill
"""
from gettext import find
from re import finditer
from sys import path
from turtle import color
import numpy as np 
import matplotlib.pyplot as plt
#rom scipy.interpolate import interp1d
import input as I
import pandas as pd
from scipy.optimize import fsolve

# parameters to define the conic shapes of the dome (all parameters defined in the paper)
c1 = I.c1
c2 = I.c2
k1 = I.k1
k2 = I.k2
h1 = I.h1
h2 = I.h2
p = I.p
er = I.er
n2 = I.n2 #dielectric refractive index
n1 = I.n1 #air refractive indeix 
wv = I.wv # wavelength in mm (defined in the paper)
k0 = I.k0 #propagation constant in free space
N = I.N
L = I.L
Array = I.Array
output_angle = I.output_angle

#const is the aribtrary variable used to center the face distribution to 0
if output_angle == 0: const = 285
elif output_angle == 20: const = 243
elif output_angle == 40: const = 440
elif output_angle == 60: const = 415
elif output_angle == 80: const = 0

m_max = 1000000000
long_r3 = h2*3
MAX_ITERATIONS = I.MAX_ITERATIONS

#=============================================================================
def f(hi, ci, ki, p): #defining the surface shapes as conics
    return hi + (ci*pow(p,2))/(1+np.sqrt(1-(1+ki)*pow(ci,2)*pow(p,2)))
    
#=============================================================================
def g(hi, ci, ki, p):
    if ci < -0.001: return 250
    else: return 350

#=============================================================================
def snell(theta_inc, n1, n2):
    arg = n1/n2 * np.sin(theta_inc)
    if abs(arg) <= 1:
        theta_ref = np.arcsin(n1 / n2 * np.sin(theta_inc))
    else:
        theta_ref = 0.
    return theta_ref
#=============================================================================


#=============================================================================
def __comp_der(f,z):
    #h = 1.E-12*z
    h = 1.E-6
    f1 = f(z-h)
    f2 = f(z+h)
    #return (f2-f1)/(2.*h)
    return np.gradient(np.array([f1, f(z), f2]), h)[1]
#=============================================================================

#=============================================================================
#def findNormal(x,ci,ki,hi):
def findNormal(x,f):
    #def F(t): return f(hi, ci, ki, t)
    def F(t): return f(t)
    m_t = __comp_der(F, x)
    if abs(m_t) == 0:
        m_n = 1.E5
    else:    
        m_n = -1./(m_t) #the slope of the normal line
    return m_n
#=============================================================================

# #=============================================================================
# def getTheta_i_max(m_n, m_n2, angle_in):
#     theta_diel2 = np.arcsin(n1/n2) #get the critical angle inside the dielectric with respect to normal 2
#     angle_btw_normals = abs(getTheta_btw(m_n, m_n2))   #get the angle between the two normals
#     theta_diel1 = - angle_btw_normals + theta_diel2 #get the critical angle inside the dielectric with respect to normal 1
#     theta_critical = np.arcsin(n2/n1*np.sin(theta_diel1))
#     theta_critical_x = getTheta_btw(0, m_n)  -  theta_critical #get the critical angle witht respect to the x axis
#     if angle_in < 0: theta_critical_x = np.pi + getTheta_btw(0, m_n) + theta_critical
#     return theta_critical_x
# #=============================================================================

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

#=============================================================================
def getUnitVector(x1,y1, x2, y2):
    vector=[x2-x1, y2-y1]
    norma = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    u = vector/norma
    return(u)
#=============================================================================

#=============================================================================
def getAmplitude(Pk, Pk1, Pk_ap, Pk_ap1, theta): #get the amplitude of the E field at the aperture plane.
    dLk = distance(Pk, Pk1)/2
    dck_ap = distance(Pk_ap, Pk_ap1)/2
    return np.sqrt(dLk/(dck_ap*np.cos(theta))), dck_ap
# =============================================================================

#=============================================================================
def getAngleBtwVectors(v1, v2):
    return np.arctan2( v1[0]*v2[1] - v1[1]*v2[0], v1[0]*v2[0] + v1[1]*v2[1] )
#=============================================================================

#=============================================================================
def getRs(theta_i, theta_t, n1, n2):
    return abs((n1*np.cos(theta_i)- n2*np.cos(theta_t))/(n1*np.cos(theta_i)+ n2*np.cos(theta_t)))**2
#=============================================================================

#=============================================================================
def getRp(theta_i, theta_t, n1, n2):
    return abs((n1*np.cos(theta_t)- n2*np.cos(theta_i))/(n1*np.cos(theta_t)+ n2*np.cos(theta_i5)))**2
#=============================================================================

#=============================================================================
def s1(x):
   return h1 + (c1*pow(x,2))/(1+np.sqrt(1-(1+k1)*pow(c1,2)*pow(x,2)))
    #return 250 + x*0
def s2(x):
   return h2  + (c2*pow(x,2))/(1+np.sqrt(1-(1+k2)*pow(c2,2)*pow(x,2)))
   #return 350

def s0(x):
    return 0   
#=============================================================================


#=============================================================================
theta_out_x2 = np.deg2rad(90-output_angle)
m_t =-1./np.tan(theta_out_x2)
x_r_max = np.cos(theta_out_x2)*long_r3 + max(Array)*np.sign(theta_out_x2)
y_r_max = abs(np.sin(theta_out_x2))*long_r3
def r3_ort(x):
        return m_t*(x - x_r_max) + y_r_max
#=============================================================================


#=============================================================================
def findIntersectionv2(fun1, Pi, vi):
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

    result_s1 = fsolve(f, [-200.0, 200.0], full_output=1)
    result_s2 = fsolve(g, [-200.0, 200.0], full_output=1)
    result_s3 = fsolve(h, [-200.0, 200.0], full_output=1)
  
    intersection = [[result_s1[0][0], s1(result_s1[0][0])], [result_s2[0][0], s2(result_s2[0][0])], [result_s3[0][0], r3_ort(result_s3[0][0])]  ]
    v = [[result_s1[0][0]-Pi[0], s1(result_s1[0][0])-Pi[1]],[result_s2[0][0]-Pi[0], s2(result_s2[0][0])-Pi[1]],[result_s3[0][0]-Pi[0], r3_ort(result_s3[0][0])-Pi[1]]]     
    j = -1
    dist = 1e5
    for i in range(0, len(v)):
        #origin = np.array([Pi[0], Pi[1]])
        #plt.quiver(*origin, *v[i], color='red')
        aux = np.dot(vi, v[i])
        #aux = distance(Pi, intersection[i])     
        if dist > aux and intersection[i] != Pi and aux > 1: 
            dist=aux
            j = i
    
  #  sk = np.append(sk, v[j])
    if j == 0:
        
        return intersection[j], s1
    elif j==1:
        return intersection[j], s2      
    elif j ==2: 
        return intersection[j], r3_ort
    else: 
        return [0,0], s0                      
#=============================================================================


def ray(vi, ri, x1, y1, iterations, i, nki, ski, Pki, ray_lengthi):
        iterations = iterations + 1
        #plt.plot(p, ri(p), color = 'red', linewidth = 0.5)
        [xi, yi], f = findIntersectionv2(ri, [x1, y1],vi)

        if [xi, yi] == [0,0]: return -1
        #Pk[i][iterations].p =  [xi, yi]
        Pki= np.append(Pki,[xi, yi])
        ray_lengthi = np.append(ray_lengthi, distance([x1, y1], [xi, yi]))

       # plt.plot(xi, yi ,'rx')
        origin = np.array([x1, y1])

        m_n = findNormal(xi, f) #find the normal of surface 1 in the intersection point 1
        v_n = np.array([1,m_n])*np.sign(m_n) #normal vector
        v_n_norm = v_n/np.sqrt(v_n[0]**2 + v_n[1]**2) #normal unit vector
        origin = np.array([xi, yi])
        theta_i = getAngleBtwVectors(v_n, vi)
        if f == s1:
            theta_r = snell(theta_i, n1, n2) #get angle_out with respect to the normal
        else: theta_r = snell(theta_i, n2, n1)   
        u = np.cos(theta_r)*v_n_norm[0] - np.sin(theta_r)*v_n_norm[1]
        v = np.sin(theta_r)*v_n_norm[0] + np.cos(theta_r)*v_n_norm[1]
        v_r = np.array([u, v])

        def r_r(x):
            return (v_r[1]/v_r[0])*(x-xi)+yi
        #plt.plot(p, r_r(p), color='green', linewidth = 0.5)
        #return nk
        if iterations < MAX_ITERATIONS:
            nki = v_n_norm
            return ray(v_r, r_r, xi, yi, iterations, i, nki, ski, Pki, ray_lengthi)   
        else: 
            ski = v_r
            return nki, ski, Pki, ray_lengthi                   
  


def directRayTracingRec(theta_i_y):
    nk = np.zeros([N,2])
    sk = np.zeros([N,2])
    row = []
    Pk = [list(row) for i in range( 0, N)]
    Ak_ap = np.zeros(N-2)
    theta_k = np.zeros(N)
    dck = np.zeros(N-2)
    phi_a = np.zeros(N)
    ray_length = [ list(row) for i in range( 0, N)]
    path_length = np.zeros(N)


    for i in range(0,len(Array)):
        theta_i_x = np.deg2rad(90-theta_i_y[i])
        x1 = Array[i]
        y1 = 0
        v1 = np.array([np.cos(theta_i_x), np.sin(theta_i_x)])
        def r1(x):
            return (v1[1]/v1[0])*(x-x1)+y1 

        #plt.plot(p, r1(p), color = 'red', linewidth = 0.5)
        iterations=0
        Pk[i]= np.append(Pk[i],[x1, y1])
        nk[i], sk[i], Pk[i], ray_length[i]  = ray(v1, r1, x1, y1, iterations, i, nk[i], sk[i], Pk[i], ray_length[i])


        d1 = ray_length[i][0]
        d2 = ray_length[i][1]   
        d3 = ray_length[i][2]
        # calculate the phase distribuiton
        phi_i = getPhaseDisrt_i(d1, d2, d3) #phase contribution due to the ray propagation
        phi_a[i] = -phi_i + const #50 is an arbitrary constant to center the phase to 0
        
        deltai = -phi_a[i]/k0
        d1 = d1 - deltai
        theta_k[i] = getAngleBtwVectors(nk[i],sk[i]) #angle between normal and pointing
        path_length[i] = d1+np.sqrt(er)*d2

        if i>1: #calculating the amplitudes
            Pstart1 = [Pk[i-2][0], Pk[i-2][1]]
            Pstart2 = [Pk[i][0], Pk[i][1]]
            Pap1 = [Pk[i-2][(MAX_ITERATIONS-1)*2], Pk[i-2][(MAX_ITERATIONS-1)*2+1]]
            Pap2 = [Pk[i][(MAX_ITERATIONS-1)*2], Pk[i][(MAX_ITERATIONS-1)*2+1]]
            Ak_ap[i-2], dck[i-2]  = getAmplitude(Pstart1, Pstart2, Pap1, Pap2, theta_k[i-2])  
            
    return Pk, Ak_ap, path_length, nk, sk, dck
    




