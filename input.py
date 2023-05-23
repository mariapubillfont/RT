# -*- coding: utf-8 -*-
"""
Created on Fri May 13 13:03:37 2022
@author: maria
input file by Marietaaaa :)


"""
import numpy as np
# parameters to define the conic shapes of the dome (all parameters defined in the paper)
c1 = -0.0021*1e3
c2 = -0.0005*1e3
k1 = -1.2
k2 = -3.9
h1 = 0.325
h2 = 0.345
D = 1.500 

p = np.linspace(-D, D, 10000)                       # number of points, represents the x-axis
e0 = 8.8541878128e-12                               #vacuum permitivitty
f = 13e9                                            #frequency            
c0 = 299792458                                      #vacuum light speed
wv = c0/f                                           # wavelength in mm (defined in the paper)
k0 = 2*np.pi/wv                                     #propagation constant in free space
L = 3*h1                                         #length of the Array (hmax = L/3) (defined in the paper)
m_max = 10000000                                    #max slope possible
const=40                                            #constant to make even phase front at aperture plane
long = 0.300                                        #how long is the final point of the rays

N = 100                                         #number of rays
Array = np.linspace (-L/2, L/2, N)                  #the starting points of the rays over the array
output_angle = 60 


tilted_angle = 20                                #in degrees from z-axis
spacing = 10                                    #ray spacing on the x-axis    
MAX_ITERATIONS = 5

#type_surface = 'flat'                               #Surface type#
#type_surface = 'conic'
type_surface = 'oblique'
#type_surface = 'tilted'


if type_surface == 'flat': h2 = 0.425
ITU_model = 1
cascade = 1
amplitude_mod = 0

matchingLayers = True
if matchingLayers:
    nSurfaces = 4 
else:
    nSurfaces = 2 

er = 2.5
mur = 1

losses = 0                                         #losses true or false
reflections = 1                                     #reflections true or false
tan_delta = 0.00066 if losses == 1 else 0           #loss coefficient in dielectric


# if reflections == 0:
#     er = np.sqrt(er)
#     mur = np.sqrt(er)


n2 = np.sqrt(er) #dielectric refractive index
n1 = 1 #air refractive indeix 
er_ML = np.sqrt(er)
nML = np.sqrt(er_ML)

n_diec = np.sqrt(er)*np.sqrt(mur)
permittivity = n_diec*n_diec*(1-1j*tan_delta)

thickness_ML1 =  wv/(np.sqrt(er_ML)*4)/1
line1_pointA = [-0.099, 0.350]
line1_pointB = [0.522, 0.074]
line2_pointA = [0.0796, 0.3722]
line2_pointB = [0.7277, 0.250]
m_line1 = (line1_pointB[1] - line1_pointA[1])/(line1_pointB[0] - line1_pointA[0])
m_line2 = (line2_pointB[1] - line2_pointA[1])/(line2_pointB[0] - line2_pointA[0])

m_tilted = np.tan(np.deg2rad(90-tilted_angle))
m_t_titled = -1./m_tilted
x_tilted0 = D
y_tilted0 = 0
thickness_titled = 0.1/np.cos(np.deg2rad(tilted_angle))

theta_out_x = np.deg2rad(90-output_angle)
theta_i_y = np.ones(N)*output_angle
y_r_max = np.sin(theta_out_x)*h2*2
x_r_max = np.cos(theta_out_x)*(D)*np.sign(output_angle)
m3 = np.tan(theta_out_x)
if m3 > m_max: m3=m_max
m_t = -1./m3
def aperture_plane(x):
    return m_t*(x - x_r_max) + y_r_max


def s0(x):
    return 0

def s1(x):
    if type_surface == 'flat' : return h1
    elif type_surface == 'conic': return h1 + (c1*pow(x,2))/(1+np.sqrt(1-(1+k1)*pow(c1,2)*pow(x,2)))
    elif type_surface == 'oblique' : return  m_line1*(x-line1_pointA[0]) + line1_pointA[1]  
    elif type_surface == 'tilted': return m_t_titled*(x-x_tilted0) + y_tilted0
    
def matchingLayer1(x):
    if type_surface == 'flat': return h1 - thickness_ML1
    elif type_surface == 'conic': return h1 - thickness_ML1 + (c1*pow(x,2))/(1+np.sqrt(1-(1+k1)*pow(c1,2)*pow(x,2)))
    elif type_surface == 'oblique' : return  m_line1*(x-line1_pointA[0]) + line1_pointA[1] - thickness_ML1/np.cos(np.arctan(m_line1))
    elif type_surface == 'tilted': return m_t_titled*(x-x_tilted0) + y_tilted0 - thickness_ML1/np.cos(np.deg2rad(tilted_angle))


def s2(x):
    if type_surface == 'flat' : return h2
    elif type_surface == 'conic': return h2  + (c2*pow(x,2))/(1+np.sqrt(1-(1+k2)*pow(c2,2)*pow(x,2)))
    elif type_surface == 'oblique' : return  m_line2*(x-line2_pointA[0]) + line2_pointA[1]
    elif type_surface == 'tilted': return m_t_titled*(x-x_tilted0) + y_tilted0 + thickness_titled 

def matchingLayer2(x):
    if type_surface == 'flat': return h2 + thickness_ML1
    elif type_surface == 'conic': return h2 + thickness_ML1 + (c2*pow(x,2))/(1+np.sqrt(1-(1+k2)*pow(c2,2)*pow(x,2)))
    elif type_surface == 'oblique' : return  m_line2*(x-line2_pointA[0]) + line2_pointA[1] + thickness_ML1/np.cos(np.arctan(m_line2))
    elif type_surface == 'tilted': return m_t_titled*(x-x_tilted0) + y_tilted0 + thickness_titled + thickness_ML1/np.cos(np.deg2rad(tilted_angle))



if type_surface == 'flat':
    surface1 = np.ones(len(p))*s1(p)
    surface2 = np.ones(len(p))*s2(p)
    MLayer1 = np.ones(len(p))*matchingLayer1(p)
    MLayer2 = np.ones(len(p))*matchingLayer2(p)
else:
    surface1 = s1(p)
    surface2 = s2(p)
    MLayer1 = matchingLayer1(p)
    MLayer2 = matchingLayer2(p)

# Find where the surface crosses z=0 and set the rest to 0
surface1 = np.where(surface1>0, surface1, 0.)
MLayer1 = np.where(MLayer1>0, MLayer1, 0.) 
surface2 = np.where(surface2>0, surface2, 0.)
MLayer2 = np.where(MLayer2>0, MLayer2, 0.)    