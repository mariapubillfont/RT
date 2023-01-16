# -*- coding: utf-8 -*-
"""
maria pubill
"""
import numpy as np 
import matplotlib.pyplot as plt
import input as I
import pandas as pd
from scipy.optimize import fsolve
import P2040_multilayer_matrix as multilayer  
import reflections_ITU_model as itu
import reverse_rayTracing

wv = I.wv # wavelength in mm (defined in the paper)
k0 = I.k0 #propagation constant in free space
N = I.N
L = I.L
Array = I.Array
output_angle = I.output_angle
p = I.p
er = I.er
mur = I.mur
e0 = I.e0
n2 = I.n_diec #dielectric refractive index
n1 = I.n1 #air refractive indeix 
er_ML = I.er_ML
nML = I.nML

type_surface = I.type_surface
thickness_ML1 = I.thickness_ML1

s1 = I.s1
s2 = I.s2
matchingLayer1 = I.matchingLayer1
matchingLayer2 = I.matchingLayer2
    

#const is the aribtrary variable used to center the face distribution to 0
if output_angle == 0: const = 0#-92
elif output_angle == 20: const = 0
elif output_angle == 40: const = 386#+36
elif output_angle == 60: const = 0
elif output_angle == 80: const = 0

m_max = 1000000000
long_r3 = I.h2*3
MAX_ITERATIONS = I.MAX_ITERATIONS



#=============================================================================
def f(hi, ci, ki, p): #defining the surface shapes as conics
    return hi + (ci*pow(p,2))/(1+np.sqrt(1-(1+ki)*pow(ci,2)*pow(p,2)))
    
def g(hi,x):
    return np.sqrt(hi**2-x**2)

#=============================================================================
def snell(theta_inc, n1, n2):
    arg = abs(n1)/abs(n2) * np.sin(theta_inc)
    if abs(arg) <= 1:
        theta_ref = np.arcsin(abs(n1) / abs(n2) * np.sin(theta_inc))
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

##=============================================================================
def getPhaseDisrt_i(d1,d2,d3,d4,d5): #get the phase distribution at the aperture plane
    return (d1+d5)*k0 + d3*k0*n2 + (d2+d4)*k0*(nML)
##=============================================================================


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

class Surface:
    def __init__(self, n1, n2, function):
        self.n1 = n1
        self.n2 = n2
        self.f = function

#=============================================================================
if I.nSurfaces == 4:
    surface1 = Surface(nML, n2, s1)
    MLayer1 = Surface(n1, nML, matchingLayer1)    
    surface2 = Surface(n2, nML, s2)
    MLayer2 = Surface(nML, n1, matchingLayer2)    
def s0(x):
    return 0   
surface0 = Surface(n1, n1, s0) 
#=============================================================================

#=============================================================================
theta_out_x2 = np.deg2rad(90-output_angle)
m_t =-1./np.tan(theta_out_x2)
x_r_max = np.cos(theta_out_x2)*long_r3 + max(Array)*np.sign(theta_out_x2)
y_r_max = abs(np.sin(theta_out_x2))*long_r3-0.600
def r3_ort(x):
        return m_t*(x - x_r_max) + y_r_max
aperture_plane = Surface(n1, n1, r3_ort)
#=============================================================================
   
def findInt(f1, f2):
    def f(xy):    
        x, y =xy 
        return np.array([y-f1(x), y-f2(x)])
    return fsolve(f, [-200.0, 200.0])
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

    result_s1 = fsolve(f, [-0.200, 0.200], full_output=1)
    result_s2 = fsolve(g, [-0.200, 0.200], full_output=1)
    result_s3 = fsolve(h, [-0.200, 0.200], full_output=1)
    result_ml1 = fsolve(ml1, [-0.200, 0.200], full_output=1)
    result_ml2 = fsolve(ml2, [-0.200, 0.200], full_output=1)
  
    intersection = [[result_ml1[0][0], matchingLayer1(result_ml1[0][0])], [result_s1[0][0], s1(result_s1[0][0])], [result_s2[0][0], s2(result_s2[0][0])], \
         [result_ml2[0][0], matchingLayer2(result_ml2[0][0])], [result_s3[0][0], r3_ort(result_s3[0][0])]]
    v = [[result_ml1[0][0]-Pi[0], matchingLayer1(result_ml1[0][0])-Pi[1]], [result_s1[0][0]-Pi[0], s1(result_s1[0][0])-Pi[1]],[result_s2[0][0]-Pi[0], s2(result_s2[0][0])-Pi[1]] , \
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
        if j == 1:
            return intersection[j], surface1
        elif j==2:
            return intersection[j], surface2      
        elif j ==4: 
            return intersection[j], aperture_plane
        elif j == 0:
            return intersection[j], MLayer1
        elif j == 3:
            return intersection[j], MLayer2         
        else: 
            return [0,0], surface0                      
#=============================================================================




def ray(vi, ri, x1, y1, iterations, i, nki, ski, Pki, ray_lengthi, all_normalsi, incident_anglei):
        iterations = iterations + 1
        #plt.plot(p, ri(p), color = 'red', linewidth = 0.5)
        [xi, yi], solution = findIntersectionv2(ri, [x1, y1],vi, 0)
        f = solution.f
        n_in = solution.n1
        n_out = solution.n2

        if [xi, yi] == [0,0]: 
            return -1
        #Pk[i][iterations].p =  [xi, yi]
        Pki= np.append(Pki,[xi, yi])
        ray_lengthi = np.append(ray_lengthi, distance([x1, y1], [xi, yi]))

        #plt.plot(xi, yi ,'rx')
        origin = np.array([x1, y1])

        m_n = findNormal(xi, f) #find the normal of surface 1 in the intersection point 1
        #all_normalsi = np.append(all_normalsi, m_n)

        v_n = np.array([1,m_n])*np.sign(m_n) #normal vector
        v_n_norm = v_n/np.sqrt(v_n[0]**2 + v_n[1]**2) #normal unit vector
        origin = np.array([xi, yi])
        theta_i = getAngleBtwVectors(v_n, vi)
        theta_t = snell(theta_i, n_in, n_out)

        u = np.cos(theta_t)*v_n_norm[0] - np.sin(theta_t)*v_n_norm[1]
        v = np.sin(theta_t)*v_n_norm[0] + np.cos(theta_t)*v_n_norm[1]
        v_t = np.array([u, v])
        #plt.quiver(*origin, *v_t, color='r')


        def r_t(x):
            return (v_t[1]/v_t[0])*(x-xi)+yi
        #plt.plot(p, r_r(p), color='green', linewidth = 0.5)
        #return nk
        if iterations < MAX_ITERATIONS:
            nki = v_n_norm
            all_normalsi = np.append(all_normalsi, m_n)
            incident_anglei = np.append(incident_anglei, theta_i)

            return ray(v_t, r_t, xi, yi, iterations, i, nki, ski, Pki, ray_lengthi, all_normalsi, incident_anglei)   
        else: 
            ski = v_t
            return nki, ski, Pki, ray_lengthi, all_normalsi, incident_anglei               
  


def directRayTracingRec(theta_i_y):
    nk = np.zeros([N,2])
    sk = np.zeros([N,2])
    row = []
    Pk = [list(row) for i in range( 0, N)]
    all_normals = [list(row) for i in range( 0, N)]
    Ak_ap = np.zeros(N-2)
    theta_k = np.zeros(N)
    incident_angle = [list(row) for i in range( 0, N)]
    dck = np.zeros(N-2)
    phi_a = np.zeros(N)
    ray_length = [ list(row) for i in range( 0, N)]
    path_length = np.zeros(N, dtype=np.complex_)
    T_coeff = np.ones(N, dtype=np.complex_)
    R_coeff= np.ones(N, dtype=np.complex_)

   

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
        nk[i], sk[i], Pk[i], ray_length[i], all_normals[i], incident_angle[i] = ray(v1, r1, x1, y1, iterations, i, nk[i], sk[i], Pk[i], ray_length[i], all_normals[i], incident_angle[i])
        d1 = ray_length[i][0]
        d2 = ray_length[i][1]   
        d3 = ray_length[i][2]
        d4 = ray_length[i][3]
        d5 = ray_length[i][4]
        theta_k[i] = getAngleBtwVectors(nk[i],sk[i]) #angle between normal and pointing


        x_in = Pk[i][2]
        y_in = Pk[i][3]
        x_ml1 = Pk[i][4]
        y_ml1  = Pk[i][5] 
        x_ml2 = Pk[i][6]
        y_ml2  = Pk[i][7]
        x_ap = Pk[i][8]
        y_ap = Pk[i][9]
        
        #v_normal = nk[i] #last normal vecto
        v_normal = np.array([1,all_normals[i][0]])*np.sign(all_normals[i][0])
        v_n_norm = v_normal/np.sqrt(v_normal[0]**2 + v_normal[1]**2) #normal unit vectorl


        if I.ITU_model:
            def r_normal(x):
                return (v_n_norm[1]/v_n_norm[0])*(x-x_in)+y_in 
            intersections = findIntersectionv2(r_normal, [x_ap, y_ap], nk[i], 1)

            layerThickness = [0, distance(intersections[0], intersections[1]), distance(intersections[2], intersections[1]), distance(intersections[2], intersections[3]), 0]
            complexPermittivity = [1, np.sqrt(er), er, np.sqrt(er), 1]
            T_coeff[i], R_coeff[i] = itu.getReflectionCoefficients_multiLayer(k0, layerThickness, 'TE', complexPermittivity, incident_angle[i])

        else:
            def r_normal(x):
                return (v_n_norm[1]/v_n_norm[0])*(x-x_in)+y_in 
            int_ML1 = findInt(r_normal, s1)  
            v_normal2 = np.array([1,all_normals[i][2]])*np.sign(all_normals[i][2])
            v_n_norm2 = v_normal2/np.sqrt(v_normal2[0]**2 + v_normal2[1]**2) #normal unit vectorl
            def r_normal2(x):
                return (v_n_norm2[1]/v_n_norm2[0])*(x-x_ml2)+y_ml2 
            int_ML2 = findInt(r_normal2, matchingLayer2)
            # #thickness1 = [distance(int_ML1, [x_in, y_in]), distance([x_ml1, y_ml1], [x_ml2, y_ml2]), distance(int_ML2, [x_ap, y_ap])]
            # plt.plot(int_ML1[0], int_ML1[1], 'bx')
            # plt.plot(x_in, y_in, 'gx')

            # plt.plot(int_ML2[0], int_ML2[1], 'black')
            # plt.plot(x_ml2, y_ml2, 'brown')
            # plt.plot(p, r_normal(p), color = 'red')
            #thickness1 = [d2, d3, d4]
            thickness1 = [distance(int_ML1, [x_in, y_in]), d3, distance(int_ML2, [x_ml2, y_ml2])]
    
            T_coeff[i] = multilayer.getReflectionCoefficients_ML(incident_angle[i], thickness1, er, I.f)

        phi_i = getPhaseDisrt_i(d1, d2, d3, d4, d5) #phase contribution due to the ray propagation
        phi_a[i] = -phi_i + const #50 is an arbitrary constant to center the phase to 0      
        deltai = -phi_a[i]/k0
        d1 = d1 - deltai

        #path_length[i] =  d1+(np.sqrt(er)-1j*np.sqrt(er)*I.tan_delta/2)*d2
        #path_length[i] =  d1+(np.sqrt(mur)*np.sqrt(er)*np.sqrt(1-1j*I.tan_delta))*d2 if I.reflections == 0 else d1+(np.sqrt(er))*d2
        path_length[i] = d1 
       
        #path_length[i] =  d1+(np.sqrt(mur)*np.sqrt(er)*np.sqrt(1-1j*I.tan_delta))*d3 + np.sqrt(er_ML)*(d2+d4) if I.reflections == 0 else d1+(np.sqrt(er))*d3 + np.sqrt(er_ML)*(d2+d4)
        #path_length[i] = d1+(np.sqrt(er))*d2

        if i>1: #calculating the amplitudes
            Pstart1 = [Pk[i-2][0], Pk[i-2][1]]
            Pstart2 = [Pk[i][0], Pk[i][1]]
            Pap1 = [Pk[i-2][(MAX_ITERATIONS-1)*2], Pk[i-2][(MAX_ITERATIONS-1)*2+1]]
            Pap2 = [Pk[i][(MAX_ITERATIONS-1)*2], Pk[i][(MAX_ITERATIONS-1)*2+1]]
            Ak_ap[i-2], dck[i-2]  = getAmplitude(Pstart1, Pstart2, Pap1, Pap2, theta_k[i-2])

    # df = pd.DataFrame(phi_a, Array)
    # df.to_excel('ph_distr_direct_' + str(output_angle) + 'deg.xlsx', sheet_name='Sheet1')


    # df_t = pd.DataFrame(abs(T_coeff), np.angle(T_coeff))
    # df_t.to_excel('t_coeff.xlsx', sheet_name='Sheet1')
    return Pk, Ak_ap, path_length, nk, sk, dck, T_coeff, R_coeff, phi_a
    