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
N =30 #number of elements of the Array
d_Array = 9 # element periodicity, in mm (defined in the paper)
d_gp = 8.4 #distance from the ground plane (defined in the paper)
Array = np.linspace (-L/2, L/2, N)
d_Array_m = Array[1] - Array[0] #measured interelement distance
theta_o_y = [80]
theta_o_x_arr = np.deg2rad([90 - x for x in theta_o_y]) #incident angle with respect to x axis
angle_in = []
m_max = 10000



# const is a constant in order to center the phase distribution to the center
if 1: #80
    const = 316/2
    
    
if 0: #40     
    const = 258.2/2  

   
if 0: #broadside
    const = 95.8

# parameters to define the conic shapes of the dome (all parameters defined in the paper)
c1 = -0.0021
c2 = -0.0005
k1 = -1.2
k2 = -3.9
h1 = 325
h2 = 345
D = 1500.
p = np.linspace(-D, D, 10000) 
er = 2.5
n2 = np.sqrt(er) #dielectric refractive index
n1 = 1. #air refractive indeix 
wv = 23. # wavelength in mm (defined in the paper)
k0 = np.pi/wv #propagation constant in free space
phi_a = []
phi_array = []
L_eff = np.zeros(len(theta_o_x_arr))
L_project = np.zeros(len(theta_o_x_arr))
M = np.zeros(len(theta_o_x_arr))
xy_min = []

#=============================================================================
def f(hi, ci, ki, p): #defining the surface shapes as conics
    return hi + (ci*pow(p,2))/(1+np.sqrt(1-(1+ki)*pow(ci,2)*pow(p,2)))
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


#==================================================
def distance(pointA, pointB):
    return (
        ((pointA[0] - pointB[0]) ** 2) +
        ((pointA[1] - pointB[1]) ** 2) 
    ) ** 0.5  # fast sqrt
#==================================================

for j in range(0,len(theta_o_x_arr)):
    theta_o_x = theta_o_x_arr[j]
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
    plt.xlabel('D (mm)')
    plt.rcParams["font.family"] = "Times New Roman"

    plt.plot(Array, np.zeros(N), '.', color='black')

    for i in range(0,len(Array)):
        ## ray 1 -> from Array to surface1 (first dome surface)
        ## ray 2 -> from surface1 to surface2 (second dome surface)
        ## ray 3 -> from surface 2 to air
            
        x1=Array[i] #points of the array
        y1=250
        #construct the line equation of ray3
        m3 = np.tan(theta_o_x)
        ray3 = m3*(p-x1)+y1
        #plt.plot(p, ray3)
       
        # find the aperture plane (perpendicular to the radiation direction)
        m_t = -1./m3
        if theta_o_y[j] >= 0:    
            x_r_max =  np.cos(theta_o_x)*h2*2 + max(Array)
        else:
            x_r_max =  np.cos(theta_o_x)*h2*2 + min(Array)
        y_r_max = abs(np.sin(theta_o_x))*h2*2 +y1
        ray3_perp =  m_t*(p - x_r_max) + y_r_max

        [xi_2,yi_2] = findIntersection(ray3, surface2, m3) #intersection between ray3 and surface2
        [xi_3,yi_3] = findIntersection(ray3, ray3_perp, m3) #intersection between ray3 and aperture plane
        plt.plot([xi_3,xi_2],[yi_3,yi_2], color='black', linewidth = 0.5) #plot the final part, from the surface 2 to the air
        
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
            print('Critical angle for element ', i+1)
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
        
        #calculate the phase distribuiton at the aperture plane
        phi_i = getPhaseDisrt_i(d1, d2, d3) #phase contribution due to the ray propagation
        
        


        #calculate the phase distribution along the central row of the illuminating array
        if abs(x0) <= max(Array):
            angle_in = np.append(angle_in, getTheta_btw(m, m_max)*180/np.pi)
            phi_a = np.append(phi_a, -phi_i + const  )
            phi_array = np.append(phi_array , x0)
        if i == 0 : [x_r_min, y_r_min] = [xi_3, yi_3]

    # L_eff[j] = distance([xy_min[0], xy_min[1]], [x_r_max, y_r_max])
    # L_project[j] = L*np.cos(theta_o_y[j]*np.pi/180)
    # M[j] = L_eff[j] / L_project[j]  
    
    # print(angle_in[0])
    # print(angle_in[len(angle_in)//2])
    # print(angle_in[len(angle_in)-1])
    
    plt.grid()
    #plt.show()
  
#plot the phase distribution
fig = plt.figure(2)
fig.set_dpi(300)
plt.plot(phi_array, angle_in)
plt.xticks([-L/2, -L/4, 0, L/4, L/2], ['-L/2', '-L/4', '0', 'L/4', 'L/2'])


fig = plt.figure(3)
fig.set_dpi(300)
ax = fig.add_subplot(111)
plt.rcParams['font.size'] = '12'
plt.plot(phi_array, phi_a)
ax.set_aspect(4, adjustable='box')

plt.yticks([-80, -40, 0, 40, 80], ['-80', '-40', '0', '40', '80'])
plt.xticks([-L/2, -L/4, 0, L/4, L/2], ['-L/2', '-L/4', '0', 'L/4', 'L/2'])
plt.ylim([-90,90])
plt.ylabel('$\phi_a$ (rad)' )
plt.xlabel('x (mm)')
plt.rcParams["font.family"] = "Times New Roman"    



plt.grid()    
# fig = plt.figure(2)
# fig.set_dpi(300)
# ax = fig.add_subplot(111)
# plt.plot(Array,phi_a[0][:], Array,phi_a[1][:], Array,phi_a[2][:])
# ax.legend([theta_o_y[0], theta_o_y[1], theta_o_y[2] ], title = '$\Theta_o$')
# ax.set_aspect(1, adjustable='box')
# plt.ylim([-90,90])
# plt.ylabel('$\phi_a$ (rad)' )
# plt.xlabel('x (mm)')
# plt.rcParams["font.family"] = "Times New Roman"    
# plt.grid()    
     
