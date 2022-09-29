# -*- coding: utf-8 -*-
"""
maria pubill
"""
from email.utils import collapse_rfc2231_value
from turtle import color
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import math
from scipy.interpolate import interp1d



# parameters to define the Array
L = 975 #length of the Array (hmax = L/3) (defined in the paper)
diagonal  = np.sqrt(pow(690,2)+pow(690,2))
N = 100 #number of elements of the Array
d_Array = 9 # element periodicity, in mm (defined in the paper)
d_gp = 8.4 #distance from the ground plane (defined in the paper)
Array = np.linspace (-L/2, L/2, N)
d_Array_m = Array[1] - Array[0] #measured interelement distance
theta_i_y = []
x = np.linspace(-L/2, L/2, 12)
xi_3_array = []
yi_3_array = []



#output angle, what you have to change.         
output_angle = 0

df = pd.read_excel('Reverse_anglesIn_' + str(output_angle) + '.xlsx', sheet_name='Sheet1')
df_np = np.array(df)
thy = df_np[:,1]
thy_array = df_np[:,0]  
const = 0
  
#plot the function of the phases
x = np.linspace(-L/2, L/2, len(thy))   
f = interp1d(x, thy, kind = 'cubic')
xnew = (np.linspace(-L/2, L/2, num=1001, endpoint=True))
theta_i_y = f(Array)
fig = plt.figure()
fig.set_dpi(300)
plt.plot(x, thy, '.')
plt.plot(xnew, f(xnew))
plt.show() 

Ak = np.ones(N)
Ak_ap = []
Pk = np.zeros([N,2])
Pk_ap = np.zeros([N,2])
dLk = []   
theta_k = []
  
    
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
k0 = 2*np.pi/wv #propagation constant in free space
phi_a = np.zeros(N)
theta_i_x_arr = np.deg2rad(90+theta_i_y)
angle_out = []
m_max = 10000


#=============================================================================
def readSurfaces():
    s1 = np.loadtxt('surface1.csv', delimiter=',')
    s2 = np.loadtxt('surface2.csv', delimiter=',')
    return s1, s2
#=============================================================================

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
    if 0:
        t = Symbol('t')
        z = hi + (ci*pow(t,2))/(1+pow(1-(1+ki)*pow(ci,2)*pow(t,2),0.5))
        dz = sym.diff(z,t) #derivate the function of the surface
        dz_value = lambdify(t,dz,'numpy')
        m_t = dz_value(x) #the slope of the tangent line m_t = f'(x), where x is the intersection of the ray with the surface
        #m_n = -1./(m_t) #the slope of the normal line    
    if 1:
        def F(t): return f(hi, ci, ki, t)
        m_t = __comp_der(F, x)
        if abs(m_t) == 0:
            m_n = 1.E5
        else:    
            m_n = -1./(m_t) #the slope of the normal line
    normal =  m_n*(p-x)+y
    #tangent =  m_t*(p-x)+y
    #plt.plot(p,normal, color = 'green', linewidth = 0.5) #if you want an auxiliar plot of the normal line
    #plt.plot(p,tangent, color='red')
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
def getAmplitude(Pk, Pk1, Pk_ap, Pk_ap1, theta): #get the amplitude of the E field at the aperture plane.
    dLk = distance(Pk, Pk1)/2
    dck_ap = distance(Pk_ap, Pk_ap1)/2
    return np.sqrt(dLk/(dck_ap*np.cos(theta)))

# =============================================================================


#defining the surfaces of the dome
if 1:
    surface1 = f(h1, c1, k1, p)
    surface1 = np.where(surface1>0, surface1, 0.)
    np.savetxt('surface1.csv', surface1, delimiter=',')
    
        
    surface2 = f(h2, c2, k2, p)
    surface2 = np.where(surface2>0, surface2, 0.)
    np.savetxt('surface2.csv', surface2, delimiter=',')
    
if 0:
    surface1 = readSurfaces()[0]
    surface2 = readSurfaces()[1]
    



    #theta_i_x = theta_i_x_arr[j]
fig = plt.figure(1)
fig.set_dpi(300)
ax = fig.add_subplot(111)
plt.plot(p, surface1, color='grey')
plt.plot(p, surface2, color='grey')
ax.set_aspect(1, adjustable='box')
ax.fill_between(p, surface1, surface2, color = 'lightgrey')
plt.ylim([0,h2*3])
plt.ylabel('z (mm)' )
plt.xlabel('x (mm)')
plt.rcParams["font.family"] = "Times New Roman"
Pk_final = []
#plt.plot(Array, np.zeros(N), 'o', color='black')

for i in range(0,len(Array)):
 ## ray 1 -> from Array to surface1 (first dome surface)
        ## ray 2 -> from surface1 to surface2 (second dome surface)
        ## ray 3 -> from surface 2 to air
        theta_i_x =  theta_i_x_arr[i]
            
        #create the line equation of the ray 1 (from the Array to surface 1)
        x1=Array[i]
        y1 = 0
        Pk[i] = [x1, y1] #save it inside an array

        m = m_max if theta_i_x == np.pi/2 else np.tan(theta_i_x)   
        m = np.tan(theta_i_x)
    
        ray1 = m*(p-x1)+y1  
        plt.plot(p, ray1, color = 'black')
        [xi,yi] = findIntersection(ray1, surface1, m)  #find the instersection between the ray1 and the surface 1 and plot ray 1   
        #plt.plot(xi, yi, 'x', color = 'red')
      # #calculate the angle_out 
        # m_n = findNormal(xi, yi, c1, k1, h1) #find the normal of surface 1 in the intersection point 1
        # theta_inc = getTheta_btw(m_n,m)
        # theta_out = snell(theta_inc, n1, n2) #get angle_out with respect to the normal
        # theta_out_x = theta_out+getTheta_btw(0,m_n)   #get angle out with respect to the x axis
       
        # #line equation that defines the ray 2 (from surface 1 to surface 2, inside the dielectric)
        # m2 = np.tan(theta_out_x)
        # ray2 = m2*(p-xi)+yi
        # #plt.plot(p,ray2, color='red')
        # [xi_2,yi_2] = findIntersection(ray2, surface2, m2) #find the instersection between the ray2 and the surface 2 and plot ray 2
        # #plt.plot(xi_2, yi_2, 'x', color = 'blue')
        # # Pk_ap[i]=[xi_2, yi_2] #points of the rays at the lens aperture (surface 2)
        # # # calculate the equation line of normal to the second surface  
        # m_n2 = findNormal(xi_2, yi_2, c2, k2, h2) #find the normal of surface 2 in the intersection point 2
    
        # # calculate the angle out
        # theta_inc2 = getTheta_btw(m_n2,m2)
        # theta_out2 = snell(theta_inc2, n2, n1) #get angle_out with respect to the normal
        # theta_out_x2 = getTheta_btw(0,m_n2) + theta_out2  #get angle out with respect to the x axis
        # if getTheta_btw(0,m_n2) < 0: #special case for negative normals
        #     theta_out_x2 = np.pi + theta_out_x2 
        # critical = getTheta_i_max(m_n, m_n2, theta_i_y[i]) #calulate the critical angle 
        # if critical > theta_i_x and theta_i_y[i] > 0 or critical < theta_i_x and theta_i_y[i] < 0 : 
        #     print('Critical angle for element ', i+1)
        #     continue
        
     
        # #line equation that defines the ray 3 (from surface 2 to air)
        # m3 = np.tan(theta_out_x2)
        # ray3 = m3*(p-xi_2)+yi_2
        # #plt.plot(p, ray3, color = 'blue')    
        # if i == 0: #case that we want an aperture plane
        #     # find the aperture plane
        #     m_t = -1./m3
        #     if theta_i_y[i] >= 0:    
        #         x_r_max =  np.cos(theta_out_x2)*h2*3 + max(Array)
        #     else:
        #         x_r_max =  np.cos(theta_out_x2)*h2*3 + min(Array)
        #     y_r_max = abs(np.sin(theta_out_x2))*h2*3 +h2
        #     ray3_perp =  m_t*(p - x_r_max) + y_r_max
        # [xi_3,yi_3] = findIntersection(ray3, ray3_perp, m3)
        # # plt.plot(p, ray3_perp)
        # x4 = xi_3
        # y4 = yi_3   

        
        # plt.plot([x1,xi],[y1,yi], color='black', linewidth = 0.5)
        # plt.plot([xi,xi_2],[yi,yi_2], color='black', linewidth = 0.5)
        # plt.plot([xi_3,xi_2],[yi_3,yi_2], color='black', linewidth = 0.5)      
        
        # #final point of the ray 3. Arbitrarly chosen, the height of this point is defined by "long"
        # # x4 = (np.cos(theta_out_x2)*long  +xi_2)
        # # y4 = abs(np.sin(theta_out_x2)*long) + yi_2
        # # Pk_final[i] = [x4, y4]       

        # angle_out = np.append(angle_out, getTheta_btw(m3, m_max)*180/np.pi)

    

    
plt.grid()
plt.show()
    
      
        
#plot the phase distribution
fig = plt.figure(2)
fig.set_dpi(300)
ax = fig.add_subplot(111)

plt.plot(Array,phi_a)
plt.title('Direct: Phase distribution for $\u03B8_o$=' + str(output_angle))
# plt.plot(Array,phi_a[0][:], Array,phi_a[1][:], Array,phi_a[2][:])
ax.set_aspect(1, adjustable='box')
plt.ylim([-90,90])
ax.set_aspect(4, adjustable='box')

plt.yticks([-80, -40, 0, 40, 80], ['-80', '-40', '0', '40', '80'])
plt.xticks([-L/2, -L/4, 0, L/4, L/2], ['-L/2', '-L/4', '0', 'L/4', 'L/2'])
# plt.ylim([-90,90])
plt.ylabel('$\phi_a$ (rad)' )
plt.xlabel('L (mm)')
plt.rcParams["font.family"] = "Times New Roman"    
plt.grid()    
     
