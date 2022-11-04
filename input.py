# -*- coding: utf-8 -*-
"""
Created on Fri May 13 13:03:37 2022
@author: maria
"""
import numpy as np


# parameters to define the conic shapes of the dome (all parameters defined in the paper)
c1 = -0.0021
c2 = -0.0005
k1 = -1.2
k2 = -3.9
h1 = 325
h2 = 425 #345

# h1 = 500
# h2 = 600

N = 100 #number of rays
D = 1500.
p = np.linspace(-D, D, 1000000) #number of points, represents the x-axis
er = 2.5 #dielectric constant of the dielectric
mur = 1
e0 = 8.8541878128e-12
n2 = np.sqrt(er) #dielectric refractive index
n1 = 1 #air refractive indeix 

f = 13e9
c0 = 299792458
wv = c0/f*1e3 # wavelength in mm (defined in the paper)
k0 = 2*np.pi/wv #propagation constant in free space
L = 3*325 #length of the Array (hmax = L/3) (defined in the paper)
Array = np.linspace (-L/2, L/2, N) #the starting points of the rays over the array
output_angle = 0 #in degrees
MAX_ITERATIONS = 3

m_max = 1000000000 #max slope possible
const=0
long = 300 #how long is the final point of the rays


losses = 1
reflections = 0

tan_delta = 0.00066 if losses == 1 else 0

if reflections == 1:
    #er = 25
    er = 2.5
    mur = 1
else:
    er = np.sqrt(2.5)
    mur = np.sqrt(2.5)    


n_diec = np.sqrt(er)*np.sqrt(mur)
permittivity = n_diec*n_diec*(1-1j*tan_delta)
alpha = np.pi*n_diec*tan_delta/wv
beta = k0 - 1j*alpha