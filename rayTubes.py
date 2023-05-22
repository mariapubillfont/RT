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

 
def getPhaseDisrt_i(ray_length, idxs, segments):                    #get the phase distribution as the wave travels
    phi_i = 0
    for j in range(0, len(idxs)):                                   #we want to include last surface (aperture)
        idx = int(idxs[j])                                          #segment index
        phi_i += ray_length[j]*segments[idx].n1*k0                  #phase 
    return phi_i

def getPathLength(rays, segments):
    #rays - all of the rays
    #segments - segments of surfaces

    path_length = np.zeros(N)
    phi_a = np.zeros(N)
    for i in range(0, len(rays)):                                   #for each ray
        ray_length = rays[i].ray_lengths                            #physical length traveled between each segment
        idxs = rays[i].idxs                                         #index for segment intersections
        phi_i = getPhaseDisrt_i(ray_length, idxs, segments)         #calculate phase shift from aperture plane to array for ray
        phi_a[i] = -phi_i                                           #initial phase on array set to get uniform phase on the output
        for j in range(0, len(idxs)-1):                             #-1 because we want to skip the last surface (aperture)
            idx = int(idxs[j])
            if  j == 0:
                path_length[i] = path_length[i] + ray_length[j]*segments[idx].n1 - phi_i/k0   
            elif reflections == 0: 
                path_length[i] = path_length[i] + ray_length[j]*segments[idx].n1

    #save phase distribution on the array
    df = pd.DataFrame(phi_a, Array)
    df.to_excel('ph_distr_direct_' + str(I.output_angle) + 'deg.xlsx', sheet_name='Sheet1')            
    return path_length

def getLastNormal(rays):
    nk =  np.zeros([N,2])                   #normal of intersection point
    sk = np.zeros([N,2])                    #poynting vector
    for i in range(0, len(rays)):
        sk[i] = rays[i].sk
        normals_aux = rays[i].normals
        nk[i] = [rays[i].normals[len(normals_aux)-4], rays[i].normals[len(normals_aux)-3]]
    return nk, sk


#=============================================================================
def calculateRayTubeAmpl(Pk, Pk1, Pk_ap, Pk_ap1, theta, theta_in):    #get the amplitude of the E field at the aperture plane
    #Pk - intersection of first ray and array
    #Pk1 - intersection of second ray and array
    #Pk_ap - intersection of first ray and aperture
    #Pk_ap1 - intersection of second ray and aperture

    dLk = distance(Pk, Pk1)/2                               #ray tube width
    dck_ap = distance(Pk_ap, Pk_ap1)/2                      #infinitesimal arc length of aperture
    if I.amplitude_mod == 1:                                
        return np.sqrt(dLk*np.cos(theta_in)/(dck_ap*np.cos(theta))), dck_ap
    else:
        return np.sqrt(dLk/(dck_ap*np.cos(theta))), dck_ap
# =============================================================================

def getAmplitude(rays, segments, theta_i):
    row = []
    Pk = [list(row) for i in range( 0, N)]                  #intersection points
    #normals = []
    Ak_ap = np.zeros(N-2)                                   #amplitude on aperture
    theta_k = np.zeros(N)
    dck = np.zeros(N-2)                                     #infinitesimal arc length of aperture    

    for i in range(0, len(rays)):                           #get all intersection points first
        Pk[i] = rays[i].Pk

    #lukas - optimize by going from 2->len(rays), maybe write with k index to follow easier
    for i in range(0, len(rays)):                                                   #for each ray
        nk = [rays[i].normals[nSurfaces*2-2], rays[i].normals[nSurfaces*2-1]]       #normal to surface
        sk = rays[i].sk                                                             #poynting vector
        theta_k[i] = getAngleBtwVectors(nk, sk)
        #theta_i[i] = np.rad2deg(rays[i].incident_angle[0]  )
        #print(theta_i[i])             
        if i > 1:                                                                   #exclude first ray, code will handle ray i-1 for each loop
            Pstart1 = [Pk[i-2][0], Pk[i-2][1]]                                      #intersection to the left of ray on array
            Pstart2 = [Pk[i][0], Pk[i][1]]                                          #intersection to the right of ray on array   
            Pap1 = [Pk[i-2][(nSurfaces)*2], Pk[i-2][(nSurfaces)*2+1]]               #intersection to the left of ray on aperture
            Pap2 = [Pk[i][(nSurfaces)*2], Pk[i][(nSurfaces)*2+1]]                   #intersection to the left of ray of aperture
            Ak_ap[i-2], dck[i-2]  = calculateRayTubeAmpl(Pstart1, Pstart2, Pap1, Pap2, theta_k[i-1], theta_i[i-1])
    return Ak_ap, dck


def getTransmissionCoef(rays, segments):
    ts_coeff = np.ones(N, dtype=np.complex_)                                                    #transmission coefficients
    ts_coeff_aggregate = np.ones(N, dtype=np.complex_)


    for i in range(0, len(rays)):
        #idxs = rays[i].idxs
        Pk = rays[i].Pk                                                                         #intersection points of ray and segment, setup as [x1 y1 x2 y2 ...]                                               
        idx = 0                                                                                 #intersection index
        intersections = np.zeros([int(len(Pk)/2)-1, 2])                                         #intersections
        thickness = []                                                                          #distance the ray travels between segments
        thickness_agrr = []
        delta = []
        incident_angle = rays[i].incident_angle
        #idxs_int = rays[i].idxs
        normals = np.zeros([len(intersections), 2])                                             #normal where there is an intersection
        #orthogonals = np.zeros([len(intersections), 2])
        orthogonal = [- rays[i].normals[1], rays[i].normals[0] ]                                #orthogonal to normal of array, counterclockwise
        last_inter = []

        #Create matrices for intersection points and corresponding surface normals
        # I could put those together and optimize it
        for j in range(0, len(Pk)-1):                                                           #-1 because we want to skip the last surface (aperture)
            if j % 2 == 0 and j > 0:                                                            #mod2 to get every intersection and j>0 to exclude array intersection
                intersections[idx] = [Pk[j], Pk[j+1]]                                           #intersection points
                idx += 1
            if j % 2 == 0 and j < len(rays[i].normals-1):                                       #mod2 to get every intersection and j>0 to exclude array intersection
                normals[idx] = [rays[i].normals[j], rays[i].normals[j+1]]                       #normals corresponding to 

        #For the ITU model the normal for the surface is extended to intersect 
        #with offseted plane to calculated travelled distance
        for j in range(0, len(intersections)-1):                                                    #for each intersection excluding aperture plane
            if j == 0:                                                                              #intersection with array, start values
                #idx_segment = int(idxs_int[j])
                v_normal = normals[j]                                                               #start normal
                #origin = intersections[j]  
                [x_0_n, y_0_n] = intersections[j]                                                   #start position
                x_end_n = x_0_n + 1*v_normal[0]                                                     #end position of ray
                y_end_n = y_0_n + 1*v_normal[1]                                                     #^
                last_inter =  [x_0_n, y_0_n]                                                        #last intersection
                #plt.plot(x_0_n, y_0_n, 'bx')
                #plt.plot(x_end_n, y_end_n, 'gx')
                #plt.quiver(*origin, *v_normal, color='red', angles='xy', scale_units='xy', scale=1)

            elif j > 0:
                thickness = np.append(thickness, distance(intersections[j], intersections[j-1]))    #distance between intersections    


                [x_0_orth, y_0_orth] = intersections[j]                                             #orhtogonal vector start
                #origin = intersections[j]
                #plt.quiver(*origin, *orthogonal, color='green', angles='xy', scale_units='xy', scale=1)
                x_end_orth = x_0_orth + 1*orthogonal[0]                                             #orthogonal vector end
                y_end_orth = y_0_orth + 1*orthogonal[1]                                             #^
               # plt.plot(x_0_orth, y_0_orth, 'o', color = 'pink')
                #plt.plot(x_end_orth, y_end_orth, 'go')

                aux = intersect_line_seg([x_0_n, y_0_n], [x_end_n, y_end_n], [x_0_orth, y_0_orth], [x_end_orth, y_end_orth])
                if aux == None:                                                                     #test if ray intersect with orthogonal vector
                    orthogonal = [ -x for x in orthogonal]                                          #rotate orthogonal vector 180 degress
                    x_end_orth = x_0_orth + 1*orthogonal[0]                                         #new orthogonal vector end test
                    y_end_orth = y_0_orth + 1*orthogonal[1]                                         #^
                    [x_int, y_int] = intersect_line_seg([x_0_n, y_0_n], [x_end_n, y_end_n], [x_0_orth, y_0_orth], [x_end_orth, y_end_orth]) #test flipped vector for intersection
                   # plt.quiver(*origin, *orthogonal, color='pink', angles='xy', scale_units='xy', scale=1)
                    #print('No intersection')
                else: 
                    [x_int, y_int] = aux
                
                #print(last_inter, [x_int, y_int])
                #
               # plt.plot(x_int, y_int, 'o', color = 'black')

                thickness_agrr = np.append(thickness_agrr, distance([x_int, y_int], last_inter))
                delta = np.append(delta, distance([x_int, y_int], [x_0_orth, y_0_orth]))
                last_inter = [x_int, y_int]    


        if I.nSurfaces == 4:
            complexPermittivity = [1, np.sqrt(er), er, np.sqrt(er), 1]
        elif I.nSurfaces == 2:
            complexPermittivity = [1, er, 1]

        #Choose transmission/reflection method
        #if I.ITU_model == 1:
        #    layerThickness =  np.pad(thickness_agrr, (1, 1), 'constant', constant_values=(0,0)) #add thickness of air
        #    ts_coeff[i] = itu.getReflectionCoefficients_multiLayer(k0, layerThickness, 'TE', complexPermittivity, incident_angle[0])
        if reflections == 1:
            ts_coeff[i] = multilayer.getReflectionCoefficients_cascade(incident_angle, thickness, er, I.f) #thickness_agg, 3rd arg
            ts_coeff_aggregate[i] = multilayer.getReflectionCoefficients_agg(incident_angle, thickness_agrr, complexPermittivity, I.f, delta)
    return ts_coeff, ts_coeff_aggregate     