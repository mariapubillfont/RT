import numpy as np
import matplotlib.pyplot as plt
from sympy import *
import input as I
import pandas as pd
import reflections_TM_model as multilayer  
import reflections_ITU_model as itu
import sys

h2 = I.h2
p = I.p
n2 = I.n_diec #dielectric refractive index
n1 = I.n1 #air refractive indeix
nML = I.nML

er = I.er
er_ML = I.er_ML 
N = I.N
L = I.L
Array = I.Array
D = I.D
m_max = I.m_max
k0 = I.k0
nSurfaces = I.nSurfaces
reflections = I.reflections

#==================================================
def print_error(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
#==================================================

#==================================================
def distance(pointA, pointB):
    return (
        ((pointA[0] - pointB[0]) ** 2) +
        ((pointA[1] - pointB[1]) ** 2)
    ) ** 0.5  # fast sqrt
#==================================================

def getAngleBtwVectors(v1, v2):
    return np.arctan2( v1[0]*v2[1] - v1[1]*v2[0], v1[0]*v2[0] + v1[1]*v2[1] )

#=============================================================================
def intersect_line_seg(p1, p2, p3, p4):
    #find the intersection between two line segments
    x1,y1 = p1
    x2,y2 = p2
    x3,y3 = p3
    x4,y4 = p4
    denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1)
    if denom == 0: # parallel
        return None
    ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom
    if ua < -1e-5 or ua > 1: # out of range
        return None
    ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom
    if ub < -1e-5 or ub > 1: # out of range
         return None
    x = x1 + ua * (x2-x1)
    y = y1 + ua * (y2-y1)
    return (x,y)
#=============================================================================


def getPhaseDisrt_i(ray_length, idxs, segments): #get the phase distribution at the aperture plane
    phi_i = 0
    for j in range(0, len(idxs)): #we want to include last surface (aperture)
        idx = int(idxs[j])
        phi_i += ray_length[j]*segments[idx].n1*k0   
    return phi_i

def getPathLength(rays, segments):
    path_length = np.zeros(N)
    phi_a = np.zeros(N)
    for i in range(0, len(rays)):
        ray_length = rays[i].ray_lengths
        idxs = rays[i].idxs
        phi_i = getPhaseDisrt_i(ray_length, idxs, segments)
        phi_a[i] = -phi_i
        for j in range(0, len(idxs)-1): #-1 because we want to skip the last surface (aperture)
            idx = int(idxs[j])
            if  j == 0:
                path_length[i] = path_length[i] + ray_length[j]*segments[idx].n1 - phi_i/k0    
            elif reflections == 0: 
                path_length[i] = path_length[i] + ray_length[j]*segments[idx].n1

    #df = pd.DataFrame(phi_a, Array)
    #df.to_excel('ph_distr_direct_' + str(I.output_angle) + 'deg.xlsx', sheet_name='Sheet1')            
    return path_length

def getLastNormal(rays):
    nk =  np.zeros([N,2])
    sk = np.zeros([N,2])
    for i in range(0, len(rays)):
        sk[i] = rays[i].sk
        normals_aux = rays[i].normals
        nk[i] = [rays[i].normals[len(normals_aux)-4], rays[i].normals[len(normals_aux)-3]]
    return nk, sk


#=============================================================================
def calculateRayTubeAmpl(Pk, Pk1, Pk_ap, Pk_ap1, theta): #get the amplitude of the E field at the aperture plane.
    dLk = distance(Pk, Pk1)/2
    dck_ap = distance(Pk_ap, Pk_ap1)/2
    return np.sqrt(dLk/(dck_ap*np.cos(theta))), dck_ap
# =============================================================================

def getAmplitude(rays, segments):
    row = []
    Pk = [list(row) for i in range( 0, N)]
    normals = []
    Ak_ap = np.zeros(N-2)
    theta_k = np.zeros(N)
    dck = np.zeros(N-2)
    for i in range(0, len(rays)):
        Pk[i] = rays[i].Pk

    for i in range(0, len(rays)):
        nk = [rays[i].normals[nSurfaces*2-2], rays[i].normals[nSurfaces*2-1]]
        sk = rays[i].sk 
        theta_k[i] = getAngleBtwVectors(nk, sk)
        if i > 1:
            Pstart1 = [Pk[i-2][0], Pk[i-2][1]]
            Pstart2 = [Pk[i][0], Pk[i][1]]
            Pap1 = [Pk[i-2][(nSurfaces)*2], Pk[i-2][(nSurfaces)*2+1]]
            Pap2 = [Pk[i][(nSurfaces)*2], Pk[i][(nSurfaces)*2+1]]
            Ak_ap[i-2], dck[i-2]  = calculateRayTubeAmpl(Pstart1, Pstart2, Pap1, Pap2, theta_k[i-2])
    return Ak_ap, dck


def getTransmissionCoef(rays, segments):
    ts_coeff = np.ones(N, dtype=np.complex_)
    row = []   
    Pk = []
    intersections = np.zeros([2, ])

    for i in range(0, len(rays)):
        idxs = rays[i].idxs
        Pk = rays[i].Pk
        idx = 0
        intersections = np.zeros([int(len(Pk)/2)-1, 2])
        thickness = []
        thickness_itu = []
        incident_angle = rays[i].incident_angle
        idxs_int = rays[i].idxs
        normals = np.zeros([len(intersections), 2])
        #orthogonals = np.zeros([len(intersections), 2])
        orthogonal = [- rays[i].normals[1], rays[i].normals[0] ]
        last_inter = []


        # I could put those together and optimize it
        for j in range(0, len(Pk)-1): #-1 because we want to skip the last surface (aperture)
            if j % 2 == 0 and j > 0:
                intersections[idx] = [Pk[j], Pk[j+1]]
                idx += 1
            if j % 2 == 0 and j < len(rays[i].normals-1):    
                normals[idx] = [rays[i].normals[j], rays[i].normals[j+1]]
                
        for j in range(0, len(intersections)-1):


            if j == 0:
                idx_segment = int(idxs_int[j])
                v_normal = normals[j]
                origin = intersections[j]
                [x_0_n, y_0_n] = intersections[j]
                x_end_n = x_0_n + 1*v_normal[0]
                y_end_n = y_0_n + 1*v_normal[1]
                last_inter =  [x_0_n, y_0_n] 
                #plt.plot(x_0_n, y_0_n, 'bx')
                #plt.plot(x_end_n, y_end_n, 'gx')
                #plt.quiver(*origin, *v_normal, color='red', angles='xy', scale_units='xy', scale=1)

            elif j > 0:
                thickness = np.append(thickness, distance(intersections[j], intersections[j-1]))


                [x_0_orth, y_0_orth] = intersections[j]
                origin = intersections[j]
                #plt.quiver(*origin, *orthogonal, color='green', angles='xy', scale_units='xy', scale=1)
                x_end_orth = x_0_orth + 1*orthogonal[0]
                y_end_orth = y_0_orth + 1*orthogonal[1]
                #plt.plot(x_0_orth, y_0_orth, 'o', color = 'pink')
                #plt.plot(x_end_orth, y_end_orth, 'go')

                aux = intersect_line_seg([x_0_n, y_0_n], [x_end_n, y_end_n], [x_0_orth, y_0_orth], [x_end_orth, y_end_orth])
                if aux == None:
                    orthogonal = [ -x for x in orthogonal]
                    x_end_orth = x_0_orth + 1*orthogonal[0]
                    y_end_orth = y_0_orth + 1*orthogonal[1] 
                    [x_int, y_int] = intersect_line_seg([x_0_n, y_0_n], [x_end_n, y_end_n], [x_0_orth, y_0_orth], [x_end_orth, y_end_orth])
                   # plt.quiver(*origin, *orthogonal, color='pink', angles='xy', scale_units='xy', scale=1)
                    #print('No intersection')
                else: 
                    [x_int, y_int] = aux
                
                #print(last_inter, [x_int, y_int])
                #plt.plot(x_int, y_int, 'o', color = 'black')

                thickness_itu = np.append(thickness_itu, distance([x_int, y_int], last_inter))
                last_inter = [x_int, y_int]    




        if I.ITU_model == 1:
            layerThickness =  np.pad(thickness_itu, (1, 1), 'constant', constant_values=(0,0))
            if I.nSurfaces == 4:
                complexPermittivity = [1, np.sqrt(er), er, np.sqrt(er), 1]
            elif I.nSurfaces == 2:
                complexPermittivity = [1, er, 1]
            else:
                print_error('Unsupported number of layers!')
            ts_coeff[i] = itu.getReflectionCoefficients_multiLayer(k0, layerThickness, 'TE', complexPermittivity, incident_angle[0])
        else:
            ts_coeff[i] = multilayer.getReflectionCoefficients_ML(incident_angle, thickness, er, I.f)
    return ts_coeff      