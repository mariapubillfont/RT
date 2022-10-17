# -*- coding: utf-8 -*-
"""
SMaria pUBILL
"""
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import openpyxl
from sympy import *
from scipy.optimize import fsolve

c1 = -0.0021
c2 = -0.0005
k1 = -1.2
k2 = -3.9
h1 = 325
h2 = 345
D = 2500.
p = np.linspace(-D, D, 20000)
er = 2.5
n2 = np.sqrt(er) #dielectric refractive index
n1 = 1. #air refractive indeix
wv = 23. # wavelength in mm (defined in the paper)
k0 = 2*np.pi/wv #propagation constant in free space
phi_a = []
phi_array = []
xy_min = []
j=0

normal1_slope = []
normal2_slope = []



# parameters to define the Array
L = h1*3 #length of the Array (hmax = L/3) (defined in the paper)
N =12 #number of elements of the Array
d_Array = 9 # element periodicity, in mm (defined in the paper)
d_gp = 8.4 #distance from the ground plane (defined in the paper)
Array = np.linspace (-L/2, L/2, N)
d_Array_m = Array[1] - Array[0] #measured interelement distance

#what you have to change
theta_o_y = 80

theta_o_x = np.deg2rad(90 -  theta_o_y) #incident angle with respect to x axis
angle_in = []
m_max = 10000
spacing= 10                                              
y1 = 0
long_r3 = h2*5 #how long is the 3rd ray (from dome surface to aperture plane)

# const is a constant used to center the phase distribution to the center (0)
if theta_o_y == 0: const = 191.5
elif theta_o_y == 20: const = 182
elif theta_o_y == 40: const = 190 
elif theta_o_y == 60: const = 0
elif theta_o_y == 80: const = 265  


# parameters to define the conic shapes of the dome (all parameters defined in the paper)
#=============================================================================
def f(hi, ci, ki, p): #defining the surface shapes as conics
    return hi + (ci*pow(p,2))/(1+np.sqrt(1-(1+ki)*pow(ci,2)*pow(p,2)))
#=============================================================================


def g(hi, ci, ki, p):
    # if ci < -0.001: return 250*np.ones(p.size)
    # else: return 350*np.ones(p.size)
    if ci < -0.001: return np.sqrt(pow(h1,2)-pow(p,2))
    else: return np.sqrt(pow(h2,2)-pow(p,2))
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
def getAngleBtwVectors(v1, v2):
    return np.arctan2( v1[0]*v2[1] - v1[1]*v2[0], v1[0]*v2[0] + v1[1]*v2[1] )
#=============================================================================

#=============================================================================
def __comp_der(f,z):
    """
    Computation of the residue of f in z
    """
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
        m_n = 1.E6
    else:
        m_n = -1./(m_t) #the slope of the normal line
    return m_n
#=============================================================================

#=============================================================================
def getTheta_i_max(  ):
    return  np.arcsin(n1/n2) #get the critical angle inside the dielectric with respect to normal 2
#=============================================================================

#=============================================================================
def getPhaseDisrt_i(d1,d2,d3): #get the phase distribution at the aperture plane
    return (d1+d3)*k0 + d2*k0*np.sqrt(er)
#=============================================================================

#=============================================================================
def crossProduct2D(u,v):
    return [u[1]*v[2]- u[2]*v[1]]
#=============================================================================

#=============================================================================
def crossProduct3D(u,v):
    return [u[2]*v[3]-u[3]*v[2], u[3]*v[1]-u[1]*v[3], u[1]*v[2]- u[2]*v[1]]
#=============================================================================


#defining the surfaces of the dome
surface1 = f(h1, c1, k1, p)

def s1(x):
   return h1 + (c1*pow(x,2))/(1+np.sqrt(1-(1+k1)*pow(c1,2)*pow(x,2)))

surface1 = np.where(surface1>0, surface1, 0.)

surface2 = f(h2, c2, k2, p)

def s2(x):
   return h2  + (c2*pow(x,2))/(1+np.sqrt(1-(1+k2)*pow(c2,2)*pow(x,2)))

def s0(x):
    return 0   

surface2 = np.where(surface2>0, surface2, 0.)
[surface_points, p_points] = getSurfacePoints(surface2, p)


#==================================================
def distance(pointA, pointB):
    return (
        ((pointA[0] - pointB[0]) ** 2) +
        ((pointA[1] - pointB[1]) ** 2)
    ) ** 0.5  # fast sqrt
#==================================================


def findIntersectionv2(fun1,fun2,x0):
    # z = np.array([y -m3*(x-x1) - y1, y - h1 + (c1*pow(x,2))/(1+np.sqrt(1-(1+k1)*pow(c1,2)*pow(x,2))) ])
    def f(xy):    
        x, y =xy 
        return np.array([y-fun1(x), y-fun2(x)])
    return fsolve(f, [-200.0, 200.0], full_output=1)    

fig = plt.figure(1)
fig.set_dpi(300)
ax = fig.add_subplot(111)
plt.rcParams['font.size'] = '12'
plt.plot(p, surface1, color='grey') 
plt.plot(p, surface2, color='grey')
# plt.plot(p, s1(p), color='grey')
# plt.plot(p, s2(p), color='grey')
ax.set_aspect(1, adjustable='box')
ax.fill_between(p, surface1, surface2, color = 'lightgrey')
#ax.fill_between(p, s1(p), s2(p), color = 'lightgrey')
plt.yticks([0, 250, 500, 750])
plt.ylim([-500,1000])
plt.xlim([-D, D])
plt.ylabel('z (mm)' )
plt.xlabel('x (mm)')
plt.rcParams["font.family"] = "Times New Roman"
# plt.plot(p_points, surface_points, '.')

y_r_max = np.sin(theta_o_x)*h2
x_r_max = np.cos(theta_o_x)*D*np.sign(theta_o_y)

m3 = np.tan(theta_o_x)
m_t = -1./m3
def r3_ort(x):
    return m_t*(x - x_r_max) + y_r_max
v3 =np.array([-np.cos(theta_o_x), -np.sin(theta_o_x)])
x = np.linspace(-D, D, 200000)



#for i in range(0,len(Array)):
for i in range(0, len(p_points)):

    ## ray 1 -> from Array to surface1 (first dome surface)
    ## ray 2 -> from surface1 to surface2 (second dome surface)
    ## ray 3 -> from surface 2 to air
  
    #points at the aperture plane, perpendicular to ray3 
    xi_3 = p_points[i]
    yi_3 = r3_ort(xi_3)
    origin3 = np.array([xi_3, yi_3])
    #plt.plot(xi_3, yi_3, 'x', color = 'pink')

    def r3(x):
        return (v3[1]/v3[0])*(x-xi_3)+r3_ort(xi_3)
  #  plt.plot(p, r3(p), color = 'blue', linewidth = 0.5)
    if m3 > m_max: m3=m_max


    #find the intersection between ray3 and surface2
    intersection2 = findIntersectionv2(s2, r3, 0.0)
    xi_2 = intersection2[0][0]
    yi_2 = s2(intersection2[0][0])
    #print(intersection2)
    #plt.plot(xi_2, yi_2, 'rx')


    m_n2 = findNormal(xi_2, c2, k2, h2) #find the slope of the normal of surface 2 in the intersection point 2
    v_n2 = -np.array([1,m_n2])*np.sign(m_n2) #normal vector
    v_n2_norm = v_n2/np.sqrt(v_n2[0]**2 + v_n2[1]**2) #normal unit vector
    origin2 = np.array([xi_2, yi_2])
    #plt.quiver(*origin, *v_n2, color='b')
    theta_out2 = getAngleBtwVectors(v_n2_norm, v3)
    theta_inc2 = snell(theta_out2, n2, n1) #get angle_out with respect to the normal

    recta_normal2 = m_n2*(p-xi_2) + yi_2 #auxilar normal line, debbugging
 #   plt.plot(p, recta_normal2, color = 'black', linewidth = 0.5)

    #we rotate the normal vector anticlockwise
    u = np.cos(theta_inc2)*v_n2_norm[0] - np.sin(theta_inc2)*v_n2_norm[1]
    v = np.sin(theta_inc2)*v_n2_norm[0] + np.cos(theta_inc2)*v_n2_norm[1]
    v2 = np.array([u, v])
    #line equation that defines the ray 2 (from surface 1 to surface 2, inside the dielectric)
    def r2(x):
        return (v2[1]/v2[0])*(x-xi_2)+yi_2
  #  plt.plot(p, r2(p), color = 'green', linewidth = 0.5)


    #theta = getAngleBtwVectors(v2, v_n2_norm)
    #print(theta*180/np.pi
    #print(theta_out2*180/np.pi, theta*180/np.pi, theta_inc2*180/np.pi)


    intersection1 = findIntersectionv2(s1, r2, 0.0)
    xi = intersection1[0][0]
    yi = s1(intersection1[0][0])
    #plt.plot(xi, yi, 'gx')


    # #calculate the angle_out
    m_n = findNormal(xi, c1, k1, h1) #find the normal of surface 1 in the intersection point 1
    v_n = -np.array([1,m_n])*np.sign(m_n) #normal vector
    v_n_norm = v_n/np.sqrt(v_n[0]**2 + v_n[1]**2) #normal unit vector
    origin = np.array([xi, yi])
    #plt.quiver(*origin, *v_n, color='b')


    theta_out = getAngleBtwVectors(v_n_norm, v2)
    theta_inc = snell(theta_out, n1, n2) #get angle_out with respect to the normal

    #we rotate the normal vector anticlockwise
    u = np.cos(theta_inc)*v_n_norm[0] - np.sin(theta_inc)*v_n_norm[1]
    v = np.sin(theta_inc)*v_n_norm[0] + np.cos(theta_inc)*v_n_norm[1]
    v1 = np.array([u, v])
    def r1(x):
        return (v1[1]/v1[0])*(x-xi)+yi

    

    intersection0 = findIntersectionv2(s0, r1, 0.0)
    x0 = intersection0[0][0]
    y0 = s0(intersection0[0][0])    

    if ( abs(theta_out) > getTheta_i_max()):
        # print('Critical angle for element ', i+1)
        continue
 

    if abs(x0) <= max(Array)+40:



        # plt.quiver(*origin3, *v3, color='r', scale = 20)
        # plt.quiver(*origin2, *v2, color='r', scale = 20)
        # plt.quiver(*origin, *v1, color='r', scale = 20)


        # calculate the distances of each ray
        d1 = distance([x0, y0],[xi, yi])
        d2 = distance([xi, yi],[xi_2, yi_2])
        d3 = distance([xi_2, yi_2],[xi_3, yi_3])

        # #calculate the phase distribuiton at the aperture plane
        phi_i = getPhaseDisrt_i(d1, d2, d3) #phase contribution due to the ray propagation    
        j = j+1

        plt.plot([x0,xi],[y0,yi], color='black', linewidth = 0.5)
        plt.plot([xi,xi_2],[yi,yi_2], color='black', linewidth = 0.5)
        plt.plot([xi_3,xi_2],[yi_3,yi_2], color='black', linewidth = 0.5)

        angle_in = np.append(angle_in, getAngleBtwVectors(-v1, [0,1])*180/np.pi)
        phi_a = np.append(phi_a, -phi_i + const  )
        phi_array = np.append(phi_array , x0)



plt.grid()
plt.show()

#save the input angles for Direct RT to excel
df = pd.DataFrame(angle_in, phi_array)
df.to_excel('Reverse_anglesIn_'+ str(theta_o_y)+'.xlsx', sheet_name='Sheet1')

#save the phase ditribution to excel
df2 = pd.DataFrame(phi_a, phi_array)
df2.to_excel('Phase_distribution_'+ str(theta_o_y) +'.xlsx', sheet_name='Sheet1')


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
plt.show()

