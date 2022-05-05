# -*- coding: utf-8 -*-
"""
SMaria pUBILL
"""
import numpy as np 
import matplotlib.pyplot as plt
#import sympy as sym
#from sympy import Symbol
#import math


# parameters to define the Array
L = 690 #length of the Array (hmax = L/3) (defined in the paper)
N =5 #number of elements of the Array
d_Array = 9 # element periodicity, in mm (defined in the paper)
d_gp = 8.4 #distance from the ground plane (defined in the paper)
Array = np.linspace (-L/2, L/2, N)
d_Array_m = Array[1] - Array[0] #measured interelement distance

#what you have to change
theta_o_y = 80

theta_o_x = np.deg2rad(90 -  theta_o_y) #incident angle with respect to x axis
angle_in = []
m_max = 10000
spacing=200



# const is a constant in order to center the phase distribution to the center
if theta_o_y == 80: #80
    const = 587
    y1=400
    
elif theta_o_y == 40: #40  
    const = 344
    y1= -400
   
elif theta_o_y == 0: #broadside
    const = 191.5
    y1=0

# parameters to define the conic shapes of the dome (all parameters defined in the paper)
c1 = -0.0021
c2 = -0.0005
k1 = -1.2
k2 = -3.9
h1 = 325
h2 = 345
D = 2000.
p = np.linspace(-D, D, 10000) 
er = 2.5
n2 = np.sqrt(er) #dielectric refractive index
n1 = 1. #air refractive indeix 
wv = 23. # wavelength in mm (defined in the paper)
k0 = 2*np.pi/wv #propagation constant in free space
phi_a = []
phi_array = []
xy_min = []
j=0

#=============================================================================
def f(hi, ci, ki, p): #defining the surface shapes as conics
    return hi + (ci*pow(p,2))/(1+np.sqrt(1-(1+ki)*pow(ci,2)*pow(p,2)))
#=============================================================================


#=============================================================================
def getSurfacePoints(s,p):
    array = []
    spoints=[]
    index = 0
    for i in range(0,len(s)-1):
        if(i%spacing==0): 
            spoints = np.append(spoints, s[i]) 
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
    #normal =  m_n*(p-x)+y
    #tangent =  m_t*(p-x)+y
    #plt.plot(p,normal, color = 'green', linewidth = 0.5) #if you want an auxiliar plot of the normal line
    #plt.plot(p,tangent, color='red')
    return m_n, m_t
#=============================================================================

#=============================================================================
def getTheta_i_max(  ):
    
    theta_diel2 = np.arcsin(n1/n2) #get the critical angle inside the dielectric with respect to normal 2
    # angle_btw_normals = abs(getTheta_btw(m_n, m_n2))   #get the angle between the two normals
    # theta_diel1 = - angle_btw_normals + theta_diel2 #get the critical angle inside the dielectric with respect to normal 1
    # theta_critical = np.arcsin(n2/n1*np.sin(theta_diel1))
    # theta_critical_x = getTheta_btw(0, m_n)  -  theta_critical #get the critical angle witht respect to the x axis
    # if theta_i_y < 0: theta_critical_x = np.pi + getTheta_btw(0, m_n) + theta_critical
    return theta_diel2
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
        if len(idx)==0: 
            print('buit')
            x=0
            y=0
            return x,y
        if len(idx) == 1:
            idx = idx[0]
            x = p[idx]
            y = f2[idx] 

        else: 
            x = p[max(idx)]
            y = f2[max(idx)]
    return x,y 
#=============================================================================

#defining the surfaces of the dome
surface1 = f(h1, c1, k1, p)
surface1 = np.where(surface1>0, surface1, 0.)
    
surface2 = f(h2, c2, k2, p)
surface2 = np.where(surface2>0, surface2, 0.)
zero_line = np.zeros(len(p))

# p_points = np.linspace(-D, D, int(10000/spacing)) 
[surface_points, p_points] = getSurfacePoints(surface2, p)

for i in range(0, len(surface_points)):
    if surface_points[i] != 0 and surface_points[i-1] == 0: surface_min = p_points[i]
    if  surface_points[i] != 0 and surface_points[i+1] == 0: surface_max = p_points[i]




#==================================================
def distance(pointA, pointB):
    return (
        ((pointA[0] - pointB[0]) ** 2) +
        ((pointA[1] - pointB[1]) ** 2) 
    ) ** 0.5  # fast sqrt
#==================================================


fig = plt.figure(1)
fig.set_dpi(300)
ax = fig.add_subplot(111)
plt.rcParams['font.size'] = '12'
plt.plot(p, surface1, color='grey')
plt.plot(p, surface2, color='grey')
ax.set_aspect(1, adjustable='box')
ax.fill_between(p, surface1, surface2, color = 'lightgrey')
plt.ylim([0,h2*3])
plt.ylabel('z (mm)' )
plt.xlabel('x (mm)')
plt.rcParams["font.family"] = "Times New Roman"
# plt.plot(p_points, surface_points, '.')


plt.plot(Array, np.zeros(N), '.', color='black')

xi_3_array = (-281.578, -135.764, -1.05011, 134.563, 280.378)
yi_3_array = (690, 690, 690, 690, 690)


#for i in range(0,len(Array)):
for i in range(0, len(p_points)):    
# for i in range(0, len(xi_3_array)):     
    
    ## ray 1 -> from Array to surface1 (first dome surface)
    ## ray 2 -> from surface1 to surface2 (second dome surface)
    ## ray 3 -> from surface 2 to air
    
        
    #construct the line equation of ray3
    x1 = p_points[i]

    # y1 = surface_points[i]
    aux = m_max*(p-x1)+y1
    # plt.plot(p_points[i], surface_points[i], 'x')
    # plt.plot(p, aux)
    # plt.plot(x1, y1, 'x')
    m3 = np.tan(theta_o_x)
    ray3 = m3*(p-x1)+y1
    # plt.plot(p, ray3, color='blue', linewidth = 0.5)
   
    # find the aperture plane (perpendicular to the radiation direction)
    m_t = -1./m3
    if theta_o_y >= 0:    
        x_r_max =  np.cos(theta_o_x)*h2*2 + surface_max
    else:
        x_r_max =  np.cos(theta_o_x)*h2*2 + surface_min
    y_r_max = abs(np.sin(theta_o_x))*h2*2 +y1
    
  
    ray3_perp =  m_t*(p - x_r_max) + y_r_max
    # plt.plot(x_r_max, y_r_max, 'x')
    # plt.plot(p,ray3_perp, color='green')
    [xi_2,yi_2] = findIntersection(ray3, surface2, m3) #intersection between ray3 and surface2
    if ( [xi_2,yi_2] == [0,0]): 
        print('**** There is not intersection between ray 3 and surface 2 ****')      
        continue
    
    # print(i)
    # plt.plot(xi_2, yi_2, 'x')
    [xi_3,yi_3] = findIntersection(ray3, ray3_perp, m3) #intersection between ray3 and aperture plane
    plt.plot([xi_3,xi_2],[yi_3,yi_2], color='black', linewidth = 0.5) #plot the final part, from the surface 2 to the air
    
    
    if 0:
        xi_3 = xi_3_array[i]
        yi_3 = yi_3_array[i]
        m3 = np.tan(theta_o_x)
        ray3 = m3*(p-xi_3)+yi_3
        [xi_2,yi_2] = findIntersection(ray3, surface2, m3)
        plt.plot([xi_3,xi_2],[yi_3,yi_2])
    
    
    # #calculate the incident angle inside the dielectric 
    m_n2 = findNormal(xi_2, yi_2, c2, k2, h2)[0] #find the normal of surface 2 in the intersection point 2
    theta_out2 = getTheta_btw(m_n2, m3)
    theta_inc2 = snell(theta_out2, n2, n1) #get angle_out with respect to the normal

    theta_inc_x2 = getTheta_btw(0,m_n2) + theta_inc2  #get angle out with respect to the x axis
    # if getTheta_btw(0,m_n2) < 0: #special case for negative normals
    #    theta_inc_x2 = np.pi - theta_inc_x2    
    #   print('normal and x if normal is negative ', theta_inc_x2*180/np.pi)
        
    #line equation that defines the ray 2 (from surface 1 to surface 2, inside the dielectric)
    m2 = np.tan(theta_inc_x2)
    ray2 = m2*(p-xi_2)+yi_2
    # plt.plot(p, ray2)
    
    
    [xi,yi] = findIntersection(ray2, surface1, m2) #find the instersection between the ray2 and the surface 2 and plot ray 2
    

    if [xi,yi] != [0,0]:
        plt.plot([xi,xi_2],[yi,yi_2], color='black', linewidth = 0.5) #plot between the surface 1 and the surface 2 (intersection 2)
    else:
        print('**** There is not intersection between ray 2 and surface 1 ****')
        import sys; sys.exit()    
    
    
    # #calculate the angle_out 
    m_n = findNormal(xi, yi, c1, k1, h1)[0] #find the normal of surface 1 in the intersection point 1
    theta_out = getTheta_btw(m_n,m2)

    theta_inc = snell(theta_out, n1, n2) #get angle_out with respect to the normal

    theta_inc_x = getTheta_btw(0,m_n) + theta_inc #get angle out with respect to the x axis

    #create the line equation of the ray 1 (from the Array to surface 1)
    m = m_max if theta_inc_x == np.pi/2 else np.tan(theta_inc_x)   
    ray1 = m*(p-xi)+yi  
    
    if ( abs(theta_out) > getTheta_i_max()):
        # print('Critical angle for element ', i+1)
        continue
    
    [x0,y0] = findIntersection(ray1, zero_line, m)  #find the instersection between the ray1 and the surface 1 and plot ray 1   
    if [x0,y0] != [0,0]:
        plt.plot([x0,xi],[y0,yi], color='black', linewidth = 0.5) #plot between the Array and the surface 1 (intersection 1)
    else:
        print('**** There is not intersection between the illuminated array and ray1 for N = ', i ,' ****')
        continue
        #import sys; sys.exit()
        
        
    
    # calculate the distances of each ray
    d1 = distance([x0, y0],[xi, yi])
    d2 = distance([xi, yi],[xi_2, yi_2])
    d3 = distance([xi_2, yi_2],[xi_3, yi_3])
    
    
    # #calculate the phase distribuiton at the aperture plane
    phi_i = getPhaseDisrt_i(d1, d2, d3) #phase contribution due to the ray propagation
    # phi_a[i] = -phi_i + const


    # calculate the phase distribution along the central row of the illuminating array
    if abs(x0) <= max(Array):
        j = j+1

        plt.plot([x0,xi],[y0,yi], color='red', linewidth = 0.5)
        plt.plot([xi,xi_2],[yi,yi_2], color='red', linewidth = 0.5)
        plt.plot([xi_3,xi_2],[yi_3,yi_2], color='red', linewidth = 0.5)
        
        angle_in = np.append(angle_in, getTheta_btw(m, m_max)*180/np.pi)
        phi_a = np.append(phi_a, -phi_i + const  )
        phi_array = np.append(phi_array , x0)
        
        
        if j == 1: 
            x_rmin = xi_3
            x_lmin = x0
            y_rmin = yi_3
            y_lmin = y0
            plt.plot(x_rmin, y_rmin,'x', color='black')
            plt.plot(x_lmin, y_lmin, 'x', color='blue')
        # if j == N:
        x_rmax = xi_3
        x_lmax = x0
        y_rmax = yi_3
        y_lmax = y0

    
     #calculation of the effective length, and magnification

    
plt.plot(x_lmax, y_lmax, 'x', color='green')
plt.plot(x_rmax, y_rmax, 'x', color='pink')
    
Leff = distance([x_rmin, y_rmin], [x_rmax, y_rmax]) #effective length at the aperture plane
Lproj = distance([x_lmin, y_lmin], [x_lmax, y_lmax]) #length of the array
M = Leff / (Lproj*np.cos(np.deg2rad(theta_o_y))) #magnification  



plt.grid()
#plt.show()
  

#plot input angles
fig = plt.figure(2)
fig.set_dpi(300)
plt.plot(phi_array, angle_in)
plt.xticks([-L/2, -L/4, 0, L/4, L/2], ['-L/2', '-L/4', '0', 'L/4', 'L/2'])


# plot the phase distribution
fig = plt.figure(3)
fig.set_dpi(300)
ax = fig.add_subplot(111)
plt.rcParams['font.size'] = '9'
plt.plot(phi_array, phi_a)
ax.set_aspect(4, adjustable='box')

plt.yticks([-80, -40, 0, 40, 80], ['-80', '-40', '0', '40', '80'])
plt.xticks([-L/2, -L/4, 0, L/4, L/2], ['-L/2', '-L/4', '0', 'L/4', 'L/2'])
plt.ylim([-90,90])
plt.title('Reverse: Phase distribution for $\u03B8_o$=' + str(theta_o_y))
plt.ylabel('$\phi_a$ (rad)' )
plt.xlabel('x (mm)')
plt.rcParams["font.family"] = "Times New Roman"    
plt.grid()    

