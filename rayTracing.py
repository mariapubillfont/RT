# -*- coding: utf-8 -*-
"""
maria pubill
"""
from gettext import find
from re import finditer
from turtle import color
import numpy as np 
import matplotlib.pyplot as plt
#rom scipy.interpolate import interp1d
import input as I
import pandas as pd
from scipy.optimize import fsolve

# import radPat



long = 600 #how long is the final point of the rays

    
# parameters to define the conic shapes of the dome (all parameters defined in the paper)
c1 = I.c1
c2 = I.c2
k1 = I.k1
k2 = I.k2
h1 = I.h1
h2 = I.h2

p = I.p
er = I.er
n2 = I.n_diec #dielectric refractive index
n1 = I.n1 #air refractive indeix 
wv = I.wv # wavelength in mm (defined in the paper)
k0 = I.k0 #propagation constant in free space
N = I.N
L = I.L
Array = I.Array
output_angle = I.output_angle
Z0 = 376.730313668

#const is the aribtrary variable used to center the face distribution to 0
if output_angle == 0: const = 285
elif output_angle == 20: const = 243
elif output_angle == 40: const = 302
elif output_angle == 60: const = 279
elif output_angle == 80: const = 0

# angle_out = []
m_max = 1000000000
long_r3 = h2*2



def smoothTriangle(data, degree):
    triangle=np.concatenate((np.arange(degree + 1), np.arange(degree)[::-1])) # up then down
    smoothed=[]

    for i in range(degree, len(data) - degree * 2):
        point=data[i:i + len(triangle)] * triangle
        smoothed.append(np.sum(point)/np.sum(triangle))
    # Handle boundaries
    smoothed=[smoothed[0]]*int(degree + degree/2) + smoothed
    while len(smoothed) < len(data):
        smoothed.append(smoothed[-1])
    return smoothed


# #=============================================================================
# def readSurfaces():
#     s1 = np.loadtxt('surface1.csv', delimiter=',')
#     s2 = np.loadtxt('surface2.csv', delimiter=',')
#     return s1, s2
# #=============================================================================

#=============================================================================
def f(hi, ci, ki, p): #defining the surface shapes as conics
    return hi + (ci*pow(p,2))/(1+np.sqrt(1-(1+ki)*pow(ci,2)*pow(p,2)))
    
#=============================================================================

def g(hi, ci, ki, p):
    if ci < -0.001: return 250*np.ones(p.size)
    else: return 350*np.ones(p.size)


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
def getTheta_btw(m1, m2):
    return np.arctan(((m2-m1)/(1+m2*m1))) #get the angle between two slopes
    #return np.arctan2(m2-m1,1+m2*m1)
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
def findNormal(x,ci,ki,hi):
    #def F(t): return f(hi, ci, ki, t)
    def F(t): return f(hi, ci, ki, t)
    m_t = __comp_der(F, x)
    if abs(m_t) == 0:
        m_n = 1.E5
    else:    
        m_n = -1./(m_t) #the slope of the normal line
    return m_n
#=============================================================================

#=============================================================================
def getTheta_i_max(m_n, m_n2, angle_in):
    theta_diel2 = np.arcsin(n1/n2) #get the critical angle inside the dielectric with respect to normal 2
    angle_btw_normals = abs(getTheta_btw(m_n, m_n2))   #get the angle between the two normals
    theta_diel1 = - angle_btw_normals + theta_diel2 #get the critical angle inside the dielectric with respect to normal 1
    theta_critical = np.arcsin(n2/n1*np.sin(theta_diel1))
    theta_critical_x = getTheta_btw(0, m_n)  -  theta_critical #get the critical angle witht respect to the x axis
    if angle_in < 0: theta_critical_x = np.pi + getTheta_btw(0, m_n) + theta_critical
    return theta_critical_x
#=============================================================================

#=============================================================================
def getPhaseDisrt_i(d1,d2,d3): #get the phase distribution at the aperture plane
    return (d1+d3)*k0 + d2*k0*np.sqrt(er)
#=============================================================================

#=============================================================================
def findIntersection(f1, f2, m):
    if 0:        
        x = 0
        y = 0
        rest = f1-f2
        for i in range(0,len(p)):
            if((rest[i] > 0 and rest[i+1] < 0) or (rest[i] < 0 and rest[i+1] > 0)):
                x = p[i]
                y = f1[i] + (f2[i]-f1[i])/2
                if m == m_max:
                    y = f2[i]
                break
    if 1:
        idx = np.argwhere(np.diff(np.sign(f1[:]-f2[:]))).flatten()
        if len(idx) == 1:
            idx = idx[0]
            x = p[idx]
            y = f2[idx]                           
    return x,y 
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


def s1(x):
   return h1 + (c1*pow(x,2))/(1+np.sqrt(1-(1+k1)*pow(c1,2)*pow(x,2)))

def s2(x):
   return h2  + (c2*pow(x,2))/(1+np.sqrt(1-(1+k2)*pow(c2,2)*pow(x,2)))

def s0(x):
    return 0   

def findIntersectionv2(fun1,fun2,x0):
    # z = np.array([y -m3*(x-x1) - y1, y - h1 + (c1*pow(x,2))/(1+np.sqrt(1-(1+k1)*pow(c1,2)*pow(x,2))) ])
    def f(xy):    
        x, y =xy 
        return np.array([y-fun1(x), y-fun2(x)])
    return fsolve(f, [-200.0, 200.0], full_output=1)  


def directRayTracing(surface1, surface2, theta_i_y, thy_array):
    #theta_i_x_arr = theta_i_y
    nk = np.zeros([N,2]) #normal of the aperture
    sk = np.zeros([N,2]) #pointying vector
    Ak = np.ones(N)
    Ak_ap = np.zeros(N-2)
    Pk = np.zeros([N,2])
    Pk_refl1 = np.zeros([N,4])
    Pk_refl2 = np.zeros([N, 4])
    Pk_intersection1 = np.zeros([N,2])
    Pk_ap = np.zeros([N,2])
    Pk_final = np.zeros([N,2])
    path_length = []
    theta_k = []
    dck = np.zeros(N-2)
    phi_a = np.zeros(N)
    angle_out = []


    for i in range(0,len(Array)):
        ## ray 1 -> from Array to surface1 (first dome surface)
        ## ray 2 -> from surface1 to surface2 (second dome surface)
        ## ray 3 -> from surface 2 to air
        theta_i_x = np.deg2rad(90-theta_i_y[i])
            
        #create the line equation of the ray 1 (from the Array to surface 1)
        x1=Array[i]
        y1 = 0
        Pk[i] = [x1, y1] #save it inside an array
        
        m = np.tan(theta_i_x) #slope of the first ray (from array to inner surface)


   
        v1 = np.array([np.cos(theta_i_x), np.sin(theta_i_x)])
        def r1(x):
            return (v1[1]/v1[0])*(x-x1)+y1 
        origin = np.array([x1, y1])
        plt.quiver(*origin, *v1, color='b',scale = 50)

        result1 = findIntersectionv2(s1, r1, 0.0)
        xi = result1[0][0]
        yi = s1(result1[0][0])
        #plt.plot(xi, yi, 'bx')

        Pk_intersection1[i] = [xi, yi]
    
        # #calculate the angle_out 
        m_n = findNormal(xi, c1, k1, h1) #find the normal of surface 1 in the intersection point 1
        v_n = np.array([1,m_n])*np.sign(m_n) #normal vector
        v_n_norm = v_n/np.sqrt(v_n[0]**2 + v_n[1]**2) #normal unit vector
        origin = np.array([xi, yi])
        theta_inc = getAngleBtwVectors(v_n, v1)
        theta_out = snell(theta_inc, n1, n2) #get angle_out with respect to the normal

        u = np.cos(theta_out)*v_n_norm[0] - np.sin(theta_out)*v_n_norm[1]
        v = np.sin(theta_out)*v_n_norm[0] + np.cos(theta_out)*v_n_norm[1]
        v2 = np.array([u, v])
        plt.quiver(*origin, *v2, color='blue', scale = 50)

        ur = -np.cos(theta_inc)*v_n_norm[0] - np.sin(theta_inc)*v_n_norm[1]
        vr = np.sin(theta_inc)*v_n_norm[0] - np.cos(theta_inc)*v_n_norm[1]   
        vr1 = np.array([ur, vr])
        #plt.quiver(*origin, *vr1, color='g')
        def r_refl1(x):
            return (vr1[1]/vr1[0])*(x-xi)+yi 
        intersection_refl1 = findIntersectionv2(s0, r_refl1, 0.0)
        xi0_r1 = intersection_refl1[0][0]
        yi0_r1 = s0(intersection_refl1[0][0])
        #plt.plot(xi0_r1, yi0_r1, 'rx')

        intersection_refl1 = findIntersectionv2(s1, r_refl1, 0.0)
        xi1_r1 = intersection_refl1[0][0]
        yi1_r1 = s1(intersection_refl1[0][0])
        #plt.plot(xi1_r1, yi1_r1, 'rx')

        Pk_refl1[i] = [xi0_r1, yi0_r1,xi1_r1, yi1_r1]

        #print(getAngleBtwVectors(v_n, v1)*180/np.pi, getAngleBtwVectors(-v_n, vr1)*180/np.pi )


        #line equation that defines the ray 2 (from surface 1 to surface 2, inside the dielectric)
        def r2(x):
            return (v2[1]/v2[0])*(x-xi)+yi
        #plt.plot(p, r2(p), color='green', linewidth = 0.5)    
        
  
        result2 = findIntersectionv2(s2, r2, 0.0)
        xi_2 = result2[0][0]
        yi_2 = s2(result2[0][0])
        Pk_ap[i]=[xi_2, yi_2] #points of the rays at the lens aperture (surface 2)

        m_n2 = findNormal(xi_2, c2, k2, h2) #find the normal of surface 2 in the intersection point 2
        v_n2 = np.array([1,m_n2])*np.sign(m_n2) #normal vector
        v_n2_norm = v_n2/np.sqrt(v_n2[0]**2 + v_n2[1]**2) #normal unit vector
        origin = np.array([xi_2, yi_2])
       # plt.quiver(*origin, *-v_n2, color='black')
        theta_inc2 = getAngleBtwVectors(v_n2_norm, v2)
        theta_out2 = snell(theta_inc2, n2, n1) #get angle_out with respect to the normal

        ur = -np.cos(theta_inc2)*v_n2_norm[0] - np.sin(theta_inc2)*v_n2_norm[1]
        vr = np.sin(theta_inc2)*v_n2_norm[0] - np.cos(theta_inc2)*v_n2_norm[1]   
        vr2 = np.array([ur, vr])
       # plt.quiver(*origin, *vr2, color='g')
        def r_refl2(x):
            return (vr2[1]/vr2[0])*(x-xi_2)+yi_2 
       # plt.plot(p, r_refl2(p), color = 'green', linewidth = 0.5)    
        intersection_refl2 = findIntersectionv2(s0, r_refl2, 0.0)
        xi0_r2 = intersection_refl2[0][0]
        yi0_r2 = s0(intersection_refl2[0][0])
        #plt.plot(xi0_r1, yi0_r1, 'rx')

        intersection_refl2 = findIntersectionv2(s2, r_refl2, 0.0)
        xi2_r2 = intersection_refl2[0][0]
        yi2_r2 = s2(intersection_refl2[0][0])
        #plt.plot(xi1_r1, yi1_r1, 'rx')

        Pk_refl2[i] = [xi0_r2, yi0_r2,xi2_r2, yi2_r2]
            
        #we rotate the normal vector anticlockwise
        u = np.cos(theta_out2)*v_n2_norm[0] - np.sin(theta_out2)*v_n2_norm[1]
        v = np.sin(theta_out2)*v_n2_norm[0] + np.cos(theta_out2)*v_n2_norm[1]
        v3 = np.array([u, v])
        #line equation that defines the ray 2 (from surface 1 to surface 2, inside the dielectric)
        def r3(x):
            return (v3[1]/v3[0])*(x-xi_2)+yi_2
        #plt.plot(p, r3(p), color = 'pink', linewidth = 0.5)
            
        
        critical = getTheta_i_max(m_n, m_n2, theta_i_y[i]) #calulate the critical angle 
        if critical > theta_i_x and theta_i_y[i] > 0 or critical < theta_i_x and theta_i_y[i] < 0 : 
            print('Critical angle for element ', i+1)
            continue
        
        theta_out_x2 = getAngleBtwVectors([1,0], v3)

        if i == 0: #case that we want an aperture plane
            # find the aperture plane
            m_t = -1./(v3[1]/v3[0])
            x_r_max = np.cos(theta_out_x2)*long_r3 + max(Array)*np.sign(theta_out_x2)
            y_r_max = abs(np.sin(theta_out_x2))*long_r3
            def r3_ort(x):
                return m_t*(x - x_r_max) + y_r_max

        result3 = findIntersectionv2(r3, r3_ort, 1.0)
        xi_3 = result3[0][0]
        yi_3 = r3(result3[0][0])
        #plt.plot(p, r3_ort(p))


        origin = np.array([xi_2, yi_2])
        plt.quiver(*origin, *v3, color='b', scale = 50)
          
        
        #final point of the ray 3. Arbitrarly chosen, the height of this point is defined by "long"
        Pk_final[i] = [xi_3, yi_3]   

        angle_out = np.append(angle_out, getAngleBtwVectors(v3, [0,1])*180/np.pi)
        #print(angle_out)
        

        # calculate the distances of each ray -> calculate the phase distribution
        d1 = distance([x1, y1],[xi, yi]) 
        d2 = distance([xi, yi],[xi_2, yi_2])     
        d3 = distance([xi_2, yi_2],[xi_3, yi_3])
        
        # deltax = (deltay-y1)/m+x1
        # plt.plot(p,np.zeros(len(p)))
        if 1:  # calculate the phase distribuiton
            phi_i = getPhaseDisrt_i(d1, d2, d3) #phase contribution due to the ray propagation
            # calculate the phase distribution along the central row of the illuminating array
            phi_a[i] = -phi_i + const #50 is an arbitrary constant to center the phase to 0
        
        deltai = -phi_a[i]/k0
        d1 = d1 - deltai

        #to calculate the amplitude at the surface 2 ---------------------------------------
        yp = h2*2 #aribitrary point in the space that fullfils the equation of the normal. We need it to calculate the unitary normal vector
        xp = (yp + m_n2*xi_2 - yi_2)/m_n2 #the x coordinates that fullfils the equation of the normal
        
        nk[i] = getUnitVector(xi_2, yi_2, xp, yp) #get the unitary vector of the normal to the surface nk
        sk[i] = getUnitVector(xi_2, yi_2, xi_3, yi_3) # poyinting vector: the direction of the ray
        theta_k = np.append(theta_k, theta_out2) #angle between normal and pointing
        path_length = np.append(path_length, d1+np.sqrt(er)*d2)
        if i>1: #calculating the amplitudes
            Ak_ap[i-2], dck[i-2]  = getAmplitude(Pk[i-2], Pk[i], Pk_ap[i-2], Pk_ap[i], theta_k[i-2])          
    

    #calculation of the effective length, and magnification
    # Effective length of the rays at the aperture -> two extreme points
    Leff_max = Pk_final[np.argmax(Pk_final[:,1])]
    Leff_min = Pk_final[np.argmin(Pk_final[:,1])]
    # Two end points of the array
    L_max = [Array[np.argmax(Array)], y1]
    L_min = [Array[np.argmin(Array)], y1]
    
    Leff = distance(Leff_max.tolist(), Leff_min.tolist()) #effective length at the aperture plane
    Lproj = distance(L_max, L_min)
    return Pk, Pk_intersection1, Pk_ap, Pk_final, sk, nk, path_length, Ak_ap, dck, theta_k , angle_out, phi_a, Leff
    





