##Direct Ray Tracing by Marietaaaa :)

import numpy as np
import matplotlib.pyplot as plt
from sympy import *
import input as I
import pandas as pd

########## input import #############
h2 = I.h2
p = I.p
n2 = I.n_diec                                   #dielectric refractive index
n1 = I.n1                                       #air refractive indeix
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

########## end input import ################

########## class definitions ##############
class LineSegment:
    def __init__(self, normal0, normal1, n1, n2, A, B, u, t, isLast, isFirst):
        self.normal0 = normal0      #normal of the starting point of the segment
        self.normal1 = normal1      #normal of the end point of the segment
        self.n1 = n1                #refractive index of medium 1
        self.n2 = n2                #refractive index of second medium
        self.A = A
        self.B = B
        self.u = u
        self.t = t
        self.isLast = isLast
        self.isFirst = isFirst                

class Ray:
    def __init__(self, Pk, sk, normals, ray_lengths, idxs, incident_angle):
        self.Pk = Pk                                #all the intersection points between the initial ray and the segments
        self.sk = sk                                #value of the last pointing vector, direction of ray coming out from the dome
        self.normals = normals                      #array with all the normals for each intersection point
        self.ray_lengths = ray_lengths              #all distances travelled by the ray, including from  the dome to the aperture plane
        self.idxs = idxs                            #indexes of all intersected segments
        self.incident_angle = incident_angle
########## end class definitions ##############


#=============================================================================
def __comp_der(f,z): #compute derivative
    # f - function of surface
    # z - x-value

    h = 1.E-6                                                   #Small step      
    f1 = f(z-h)                                                 #Just to the left 
    f2 = f(z+h)                                                 #Just to the right    
    return np.gradient(np.array([f1, f(z), f2]), h)[1]          #Calculate the gradient and take out the value corresponding to f(z), returns [3,] array
#=============================================================================

#=============================================================================
def getNormal(x,f):
    # x - x-value
    # f - function for surface
  
    m_t = __comp_der(f, x)                  #tangential slope
    if abs(m_t) == 0:                       
        m_n = 1.E5                          #normal slope is infinite, approx 10^5
    else:    
        m_n = -1./(m_t)                     #the slope of the normal line
    return m_n
#=============================================================================

#=============================================================================
def discretize_function(f, n_segments, n1, n2, isLast, isFirst):
    # f - the function we want to discretize in n segments.
    # n_segments - number of segments to discretize
    # n1 - refractive index of the first medium
    # n2 - refractive index of second medium
    # isLast - Bool, indicates if the segment is the furthest out
    # isFirst - Bool, indicates if segment is the array

    surface = []                                                                            #initiate surface vactor     
    discrete_x = np.linspace(-D, D, n_segments+1)                                           #start and stop x-values for segments
    for i in range(0, n_segments):
        u = [discrete_x[i+1] - discrete_x[i], f(discrete_x[i+1]) - f(discrete_x[i])]        #maps out the straight line that is the segment

        m_n0 = getNormal(discrete_x[i], f)
        v_n0 = np.array([1,m_n0])*np.sign(m_n0)                                             #normal vector of first point
        v_unitn0 = v_n0/np.sqrt(v_n0[0]**2 + v_n0[1]**2)                                    #normal unit vector of first point

        m_n1 = getNormal(discrete_x[i+1], f)
        v_n1 = np.array([1,m_n1])*np.sign(m_n1)                                             #normal vector of second point
        v_unitn1 = v_n1/np.sqrt(v_n1[0]**2 + v_n1[1]**2)                                    #normal unit vector of second point

        surface = np.append(surface, LineSegment(v_unitn0, v_unitn1, n1, n2, [discrete_x[i], f(discrete_x[i])], [discrete_x[i+1], f(discrete_x[i+1])], u, 1, isLast, isFirst))
    return surface
#=============================================================================

#=============================================================================
#calculate points where the ray intersects the surface
def getSurfacePoints(s,p):
    # s - vector of points corresponding to surface s(p)
    # p - vector of x-points 
    array = []
    spoints=[]
    s_aux = np.ones(len(p))*s if type(s) == float else s
    for i in range(0,len(s_aux)-1):
        if(i%I.spacing==0):                                         #find where the rays start
            spoints = np.append(spoints, s_aux[i])
            array = np.append(array, p[i])

    return spoints, array
#=============================================================================


#==================================================
def distance(pointA, pointB):
    return (
        ((pointA[0] - pointB[0]) ** 2) +
        ((pointA[1] - pointB[1]) ** 2)
    ) ** 0.5  # fast sqrt
#==================================================


#=============================================================================
def getAngleBtwVectors(v1, v2):
    return np.arctan2( v1[0]*v2[1] - v1[1]*v2[0], v1[0]*v2[0] + v1[1]*v2[1] )
#=============================================================================

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
#lukas - i does nothing
def intersect_line_seg(segm, i, p3, p4):
    #segm - segment of surface to be checked
    #i - supposedly the index of segment but it does nothing
    #p3 - ray start coordinates
    #p4 - ray end values

    #find the intersection between two line segments
    x1,y1 = segm.A                                              #Start values for segment
    x2,y2 = segm.B                                              #End values for segment    
    x3,y3 = p3                                                  #ray start values
    x4,y4 = p4                                                  #ray end values
    denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1)                   #tilt difference
    if denom == 0: # parallel
        return None
    ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom
    if ua < -1e-5 or ua > 1: # out of range
        return None
    ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom
    if ub < -1e-5 or ub > 1: # out of range
         return None
    x = x1 + ua * (x2-x1)                                       #Intersection point x
    y = y1 + ua * (y2-y1)                                       #Intersection point y
    return (x,y), i
#=============================================================================

#=============================================================================
def get_intersection(segments, rayi_A, rayi_B, direction):
    #segments - surface segments
    #rayi_A - start point for ray
    #rayi_B - end point for ray
    #direction - direct or reverse ray tracing

    Ray = [rayi_B[0] - rayi_A[0], rayi_B[1] - rayi_A[1]]            #line that ray would shoot
    dist = 1e5                                                      #starting value for maximum distance travelled, very big 
    for i in range(0, len(segments)):                               #check each segment
        #plt.quiver(*segments[i].A, *segments[i].u, color = 'red', angles='xy', scale_units='xy', scale=1)

        aux = intersect_line_seg(segments[i], i, rayi_A, rayi_B)
        #if there is an intersection, does not calculate travelled length but which segment got the closest intersection
        if aux != None:
            intersection = aux[0]                                                           #intersection between ray and segment
            Vector = [intersection[0]-rayi_A[0], intersection[1] - rayi_A[1]]               #vector from start to next intersection point
            prod =  np.dot(Ray, Vector)                                                     #auxillary distance for comparison between intersections 
            if dist > prod and (prod >= 1e-5 or ((segments[i].isFirst and direction == 'reverse')\
                 or (segments[i].isLast and direction == 'direct') )):# or segments[i].isLast):
                dist = prod
                next_idx = i
                [xi, yi] = intersection
    if dist == 1e5:
        return None                    
    else:
        return [xi, yi], next_idx    
#=======================================121378======================================



def findThickness(segments, pointA, pointB):
    thickness = []
    inter = []
    Ray = [pointB[0] - pointA[0], pointB[1] - pointA[1]]
    idx = 0
    for i in range(0, len(segments)):
        aux = intersect_line_seg(segments[i], i, pointA, pointB)
        if aux != None:
            intersection = aux[0]
            Vector = [intersection[0]-pointA[0], intersection[1] - pointB[1]]
            prod =  np.dot(Ray, Vector)
            if (not segments[i].isLast and not segments[i].isFirst) and abs(prod) >= 1e-5:
                inter = np.append(inter, aux[0])
                plt.plot(intersection[0], intersection[1], 'gx')
                if idx > 0:
                    thickness=np.append(thickness, distance([inter[(idx-1)*2], inter[(idx-1)*2 +1 ]], [inter[(idx)*2], inter[(idx)*2 +1 ]]))
                idx += 1
    return inter, thickness    




#=============================================================================
def findInterNormal(N0, N1, xi, yi, A, B ):
    #N0 - start of segment normal
    #N1 - end of segment normal
    #xi - intersect x point
    #yi - intersect y point
    #A - segment start point
    #B - segment end point

    #https://gamedev.stackexchange.com/questions/18615/how-do-i-linearly-interpolate-between-two-vectors
    # Interpolate the normals for the intersection point
    xa = A[0]
    ya = A[1]
    xb = B[0]
    yb = B[1]
    distance_a = np.sqrt((xi-xa)**2 + (yi-ya)**2)
    distance_b = np.sqrt((xb-xi)**2 + (yb-yi)**2)
    distance_ab = np.sqrt((xb-xa)**2 + (yb-ya)**2)
    da = distance_a/distance_ab
    db = distance_b/distance_ab
    alpha = da
    interNormal = alpha*N1 + (1-alpha)*N0
    return interNormal/np.linalg.norm(interNormal)
#=============================================================================



#=============================================================================
def ray(segments, rayi_A, rayi_B, Pki, direction, ski, ray_lengthi, all_normalsi, incident_anglei, idx_segmentsi):
    #segments - list of surface segments
    #rayi_A - start values for ray
    #rayi_B - end values for ray 
    #Pki - surface intersection points for ray
    #direction - direct or reverse ray tracing
    #ski - poynting vector
    #ray_lengthi
    #all_normalsi
    #incident_anglei
    #idx_segmentsi

    #problem with direct and last surface

    # Find intersection of ray and segment
    aux = get_intersection(segments, rayi_A, rayi_B, direction)
    if aux == None: 
        #print("No intersection found")
        return None

    [x_int, y_int] = aux[0]                                                         #Intersection points
    idx_int = aux[1]                                                                #index of intersected segment    
    Pki = np.append(Pki, [x_int, y_int])                                            #append intersection point    
    Ray_in = [rayi_B[0] - rayi_A[0], rayi_B[1] - rayi_A[1]]                         #input ray definition

    Normal0 = segments[idx_int].normal0                                             #normal at segment start
    Normal1 = segments[idx_int].normal1                                             #normal at segment end
    # Interpolate to find normal at intersection point
    Normal = findInterNormal(Normal0, Normal1, x_int, y_int, segments[idx_int].A, segments[idx_int].B)  

    ray_lengthi = np.append(ray_lengthi, distance(rayi_A,[x_int, y_int]))
    idx_segmentsi = np.append(idx_segmentsi, idx_int)
    all_normalsi = np.append(all_normalsi, Normal)

    theta_i = getAngleBtwVectors(Normal, Ray_in)                                    #incoming angle of ray to the surface
    if np.dot(Normal, Ray_in) > 0:                                                  #forwdard direction, normal and ray same direction. From epsilon1 to epsilon2
        theta_t = snell(theta_i, segments[idx_int].n1, segments[idx_int].n2)
    else:                                                                           #reverse direction, normal and ray opposite directions    
        theta_t = np.pi - snell(theta_i, segments[idx_int].n2, segments[idx_int].n1)
        if (abs(theta_t == np.pi) and segments[idx_int].n2 != 1 and I.type_surface != 'flat'): #critical angle, when the surface is not flat
            #print("critical angle")
            return None  
    u = np.cos(theta_t)*Normal[0] - np.sin(theta_t)*Normal[1]                                          
    v = np.sin(theta_t)*Normal[0] + np.cos(theta_t)*Normal[1]
    Ray_t = np.array([u, v])                                                        #new ray direction
    x_end = x_int + 3*Ray_t[0]                                                      #new ray end point
    y_end = y_int + 3*Ray_t[1]                                                      #^

    #check if the ray is at the last or first surface or recursive call ray with new start and end points for ray
    if segments[idx_int].isFirst and direction == 'reverse' or segments[idx_int].isLast and direction == 'direct':
        ski = Ray_t
        return Pki, ski, ray_lengthi, all_normalsi, incident_anglei, idx_segmentsi
    else:
        #nki = Normal
        incident_anglei = np.append(incident_anglei, theta_i)
        return ray(segments, [x_int, y_int], [x_end, y_end], Pki, direction, ski, ray_lengthi, all_normalsi, incident_anglei, idx_segmentsi)
#=============================================================================    


#=============================================================================  
def directRayTracing_segments(theta_i_y, segments):
    # theta_i_y - starting angle
    #segments - surface segments

    t = 3
    row = []
    Pk = [list(row) for i in range( 0, N)]
    ray_length = [ list(row) for i in range( 0, N)]
    incident_angle = [list(row) for i in range( 0, N)]
    all_normals = [list(row) for i in range( 0, N)]
    sk = np.zeros([N,2])
    idx_segments = [list(row) for i in range( 0, N)]
    rays = []

    for i in range(0, len(Array)):                                                              #for each ray
        xA = Array[i]                                                                           #start point for ray
        yA = 0                                                                                  #^    
        xB = xA + np.sin(theta_i_y[i])*t                                                        #end point for ray
        yB = yA + np.cos(theta_i_y[i])*t                                                        #^

        Pk[i] = np.append(Pk[i], [xA, yA])                                                      #append atart intersection point

        Pk[i], sk[i], ray_length[i],all_normals[i], incident_angle[i], idx_segments[i] = \
            ray(segments, [xA, yA], [xB, yB], Pk[i], 'direct', sk[i], ray_length[i], all_normals[i], incident_angle[i], idx_segments[i])
        rays = np.append(rays, Ray(Pk[i], sk[i], all_normals[i], ray_length[i], idx_segments[i], incident_angle[i]))
    return rays
#=============================================================================  


#=============================================================================  
def reverseRayTracing_segments(theta_out_x, segments):
    # theta_out_x - angle of ray defined from x-axis
    # segments - the segments of all surfaces

    #y is used instead of z
    theta_out_y = np.pi/2 - np.deg2rad(theta_out_x)                             #angle defined from z-axis
    [y_aperture, x_aperture] = getSurfacePoints(I.aperture_plane(p), p)         #calculate starting points for the rays on the aperture
    row = []
    Pk = [list(row) for i in range( 0, len(x_aperture))]                        #list of list to handle each rays intersection points
    angle_in = []
    angle_position = []
    ray_length = [ list(row) for i in range( 0, len(x_aperture))]
    incident_angle = [list(row) for i in range( 0, len(x_aperture))]
    all_normals = [list(row) for i in range( 0, len(x_aperture))]
    idx_segments = [list(row) for i in range( 0, len(x_aperture))]
    sk = np.zeros([len(x_aperture),2])
    plot = True

    # FOR PLOTTING PRUPOSES ONLY
    if plot:            
        fig = plt.figure()
        fig.set_dpi(700)
        ax = fig.add_subplot(111)
        font = {'family' : 'Times New Roman',
                'weight' : 'normal',
                'size'   : 12}
        plt.rc('font', **font)
        plt.ylim([0, h2*2])
        plt.xlim([-D, D])
        plt.ylabel('z (m)')
        plt.xlabel('x (m)')
        plt.plot(p, surface1_arr, color='gray', linewidth = 1)
        plt.plot(p, surface2_arr, color='gray', linewidth = 1)
        ax.fill_between(p, surface2_arr, surface1_arr, color = '#c3c5e2')
        ax.set_aspect(1, adjustable='box')


    for i in range(0, len(x_aperture)):         #For each ray
        xA = x_aperture[i]                      #x starting point for ray 
        yA = y_aperture[i]                      #z starting point for ray 
        t = -yA/np.sin(theta_out_y)             
        xB = xA + np.cos(theta_out_y)*t         #x value for line where ray would intersect y=0 if no other surface
        yB = 0                                  #y value for line where ray would intersect y=0 if no other surface
        Pk[i] = np.append(Pk[i], [xA, yA])      #Starting intersection with aperture plane
         
        #Calculate the rest of the intersections points for the ray 
        result = ray(segments, [xA, yA], [xB, yB], Pk[i], 'reverse', sk[i], ray_length[i], all_normals[i], incident_angle[i], idx_segments[i])
        if result != None:
            Pk[i] = result[0]                   #Extract list of intersection points
            x0 = Pk[i][len(Pk[i])-2]            #points intersecting the array
            y0 = Pk[i][len(Pk[i])-1]            #^
            x1 =  Pk[i][len(Pk[i])-4]           #points from second to last intersection
            y1 =  Pk[i][len(Pk[i])-3]           #^
            lastRay = [x1-x0, y1-y0]            #ray that ends up at array

            # If the ray hits the array
            if abs(x0) <= max(Array) + max(Array)*0.1 and y0 >= 0:                          #add a 10% of tolerance
                angle_in = np.append(angle_in, getAngleBtwVectors(lastRay, [0,1]))          #angle between array normal and ray
                angle_position = np.append(angle_position, x0)                              #where ray hits array
                Pk_np = np.array(Pk[i])
                if plot:
                    for m in range(0,int(len(Pk_np)/2) - 1):
                       plt.plot([Pk_np[m*2], Pk_np[m*2+2]], [Pk_np[m*2+1], Pk_np[m*2+3]], color='black', linewidth = 0.5)
    if plot:
        plt.show()
    return angle_in, angle_position
#=============================================================================  



