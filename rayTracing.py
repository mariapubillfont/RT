# -*- coding: utf-8 -*-
"""
maria pubill
"""
import numpy as np 
import numba
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import input as I
# import radPat



const=0
long = 300 #how long is the final point of the rays

    
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
# angle_out = []
m_max = 1000000000



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
#=============================================================================

#=============================================================================
def __comp_der(f,z):
    """
    Computation of the residue of f in z
    """
    h = 1.E-12*z
    f1 = f(z-h)
    f2 = f(z+h)
    return (f2-f1)/(2.*h)
#=============================================================================

#=============================================================================
def findNormal(x,y,ci,ki,hi):
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
# def getRadPattern(): #calculate the radiattion pattern for
#     dLk = distance(Pk, Pk1)/2
#     dck_ap = distance(Pk_ap, Pk_ap1)/2
#     return np.sqrt(dLk/(dck_ap*np.cos(theta)))
# =============================================================================



def directRayTracing(surface1, surface2, theta_i_y):
    Array = np.linspace (-L/2, L/2, N)
    theta_i_x_arr = np.deg2rad(90+theta_i_y)
    nk = np.zeros([N,2]) #normal of the aperture
    sk = np.zeros([N,2]) #pointying vector
    Ak = np.ones(N)
    Ak_ap = np.zeros(N-2)
    Pk = np.zeros([N,2])
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
        theta_i_x = theta_i_x_arr[i]
            
        #create the line equation of the ray 1 (from the Array to surface 1)
        x1=Array[i]
        # y1=-250
        y1 = -400
        Pk[i] = [x1, y1] #save it inside an array
        
        
        
        m = m_max if theta_i_x == np.pi/2 else np.tan(theta_i_x)   
        m = np.tan(theta_i_x)
    
        ray1 = m*(p-x1)+y1   
        [xi,yi] = findIntersection(ray1, surface1, m)  #find the instersection between the ray1 and the surface 1 and plot ray 1   
        # xi = x1 if theta_i_x == np.pi/2 else xi
        # plt.plot(p,ray1,'red')
        Pk_intersection1[i] = [xi, yi]
        # print(m*(-345-x1)+y1)
    
        # #calculate the angle_out 
        m_n = findNormal(xi, yi, c1, k1, h1) #find the normal of surface 1 in the intersection point 1
        theta_inc = getTheta_btw(m_n,m)
        theta_out = snell(theta_inc, n1, n2) #get angle_out with respect to the normal
        theta_out_x = theta_out+getTheta_btw(0,m_n)   #get angle out with respect to the x axis
       
        #line equation that defines the ray 2 (from surface 1 to surface 2, inside the dielectric)
        m2 = np.tan(theta_out_x)
        ray2 = m2*(p-xi)+yi
        [xi_2,yi_2] = findIntersection(ray2, surface2, m2) #find the instersection between the ray2 and the surface 2 and plot ray 2
        Pk_ap[i]=[xi_2, yi_2] #points of the rays at the lens aperture (surface 2)
        # plt.plot(xi_2, yi_2, 'x')
        # plt.plot(p,ray2, 'green')
        # calculate the equation line of normal to the second surface  
        m_n2 = findNormal(xi_2, yi_2, c2, k2, h2) #find the normal of surface 2 in the intersection point 2
    
        # calculate the angle out
        theta_inc2 = getTheta_btw(m_n2,m2)
        theta_out2 = snell(theta_inc2, n2, n1) #get angle_out with respect to the normal
        theta_out_x2 = getTheta_btw(0,m_n2) + theta_out2  #get angle out with respect to the x axis
        if getTheta_btw(0,m_n2) < 0: #special case for negative normals
            theta_out_x2 = np.pi + theta_out_x2 
        critical = getTheta_i_max(m_n, m_n2, theta_i_y[i]) #calulate the critical angle 
        if critical > theta_i_x and theta_i_y[i] > 0 or critical < theta_i_x and theta_i_y[i] < 0 : 
            print('Critical angle for element ', i+1)
            continue
        
     
        #line equation that defines the ray 3 (from surface 2 to air)
        m3 = np.tan(theta_out_x2)
        ray3 = m3*(p-xi_2)+yi_2
        # plt.plot(p, ray3, 'blue')
    
        if 0: #case that we want an aperture plane
            # find the aperture plane
            m_t = -1./m3
            if theta_i_y[i] >= 0:    
                x_r_max =  np.cos(theta_out_x2)*h2*2 + max(Array)
            else:
                x_r_max =  np.cos(theta_out_x2)*h2*2 + min(Array)
            y_r_max = abs(np.sin(theta_out_x2))*h2*2 +y1
            ray3_perp =  m_t*(p - x_r_max) + y_r_max
            [xi_3,yi_3] = findIntersection(ray3, ray3_perp, m3)
    
        
        #final point of the ray 3. Arbitrarly chosen, the height of this point is defined by "long"
        x4 = (np.cos(theta_out_x2)*long  +xi_2)
        y4 = abs(np.sin(theta_out_x2)*long) + yi_2
        Pk_final[i] = [x4, y4]       
        
        angle_out = np.append(angle_out, getTheta_btw(m3, m_max)*180/np.pi)
        
        #calculation of the effective length, and magnification
        if 0:
            if i == 0: 
                x_rmin = x4
                x_lmin = x1
                y_rmin = y4
                y_lmin = y1
    
            if i == N-1:
                x_rmax = x4
                x_lmax = x1
                y_rmax = y4
                y_lmax = y1
                Leff = distance([x_rmin, y_rmin], [x_rmax, y_rmax]) #effective length at the aperture plane
                Lproj = distance([x_lmin, y_lmin], [x_lmax, y_lmax]) #length of the array
                # M = Leff / (Lproj*np.cos(np.deg2rad(output_angle))) #magnification  
        
        
        deltai = (L - i*(L/N))*np.cos(theta_i_x)
        # print(theta_i_x*180/np.pi)
        print(deltai)
        # calculate the distances of each ray -> calculate the phase distribution
        d1 = distance([x1, y1],[xi, yi]) - deltai
        # print(d1)
        d2 = distance([xi, yi],[xi_2, yi_2])     
        d3 = distance([xi_2, yi_2],[x4, y4])
        
        
        if 0:  # calculate the phase distribuiton
            phi_i = getPhaseDisrt_i(d1, d2, d3) #phase contribution due to the ray propagation
            # calculate the phase distribution along the central row of the illuminating array
            phi_a[i] = -phi_i + const #50 is an arbitrary constant to center the phase to 0
        
        
        #to calculate the amplitude at the surface 2 ---------------------------------------
        yp = h2*2 #aribitrary point in the space that fullfils the equation of the normal. We need it to calculate the unitary normal vector
        xp = (yp + m_n2*xi_2 - yi_2)/m_n2 #the x coordinates that fullfils the equation of the normal
        
        nk[i] = getUnitVector(xi_2, yi_2, xp, yp) #get the unitary vector of the normal to the surface nk
        sk[i] = getUnitVector(xi_2, yi_2, x4, y4) # poyinting vector: the direction of the ray
        theta_k = np.append(theta_k, getTheta_btw(m_n2, m3)) #angle between normal and pointing
        path_length = np.append(path_length, d1+np.sqrt(er)*d2)
        if i>1: #calculating the amplitudes
            
                 
            Ak_ap[i-2], dck[i-2]  = getAmplitude(Pk[i-2], Pk[i], Pk_ap[i-2], Pk_ap[i], theta_k[i-2])
            

            
    # f = interp1d(Array[1:N-1], Ak_ap, kind = 'cubic')        
    # Ak_new = smoothTriangle(Ak_ap, 2)             
    # fig = plt.figure(7)
    # fig.set_dpi(300)
    # plt.plot(Array[1:N-1], smoothTriangle(Ak_ap, 2), color='red')
    # plt.ylabel('diferencia a larray')
    # plt.xlabel('Array x(mm)')
    # plt.grid()   
    # plt.show()


  
    # Ak_ap = Ak_ap*np.cos(theta_i_x)
    return Pk, Pk_intersection1, Pk_ap, Pk_final, sk, nk, path_length, Ak_ap, dck, theta_k , angle_out 
    





