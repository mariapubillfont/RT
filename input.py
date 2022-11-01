# -*- coding: utf-8 -*-
"""
Created on Fri May 13 13:03:37 2022

@author: maria
"""
from math import perm
import numpy as np


const=0
long = 300 #how long is the final point of the rays

    
# parameters to define the conic shapes of the dome (all parameters defined in the paper)

c1 = -0.0021
c2 = -0.0005
k1 = -1.2

k2 = -3.9
h1 = 325

#h2 = 345 #345
h2 = 425

N = 2
D = 1500.
p = np.linspace(-D, D, 1000000) 
er = 25
e0 = 8.8541878128e-12

#n2 = np.sqrt(er) 
n_diec = 5
# #dielectric refractive index
n1 = 1 #air refractive indeix 
wv = 23.0769 # wavelength in mm (defined in the paper)
k0 = 2*np.pi/wv #propagation constant in free space
L = 3*h1 #length of the Array (hmax = L/3) (defined in the paper)
Array = np.linspace (-L/2, L/2, N)
#Array = (-223.2723272,99.759976,345.2845285)
#Array = (-L/2, L/2)
output_angle = 40
MAX_ITERATIONS = 3
#tan_delta = 0.00066
tan_delta = 0.00066
permittivity = n_diec*n_diec*(1-1j*tan_delta)
#n2 = permittivity

alpha = np.pi*n_diec*tan_delta/wv
beta = k0 - 1j*alpha

# angle_out = []
m_max = 1000000000

