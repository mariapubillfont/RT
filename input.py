# -*- coding: utf-8 -*-
"""
Created on Fri May 13 13:03:37 2022

@author: maria
"""
import numpy as np


const=0
long = 300 #how long is the final point of the rays

    
# parameters to define the conic shapes of the dome (all parameters defined in the paper)

c1 = -0.0021
c2 = -0.0005
k1 = -1.2

k2 = -3.9
h1 = 325

h2 = 345 #345


N = 100
D = 1500.
p = np.linspace(-D, D, 1000000) 
er = 2.5
n2 = np.sqrt(er) #dielectric refractive index
n1 = 1 #air refractive indeix 
wv = 23. # wavelength in mm (defined in the paper)
k0 = 2*np.pi/wv #propagation constant in free space
L = 3*h1 #length of the Array (hmax = L/3) (defined in the paper)
Array = np.linspace (-L/2, L/2, N)


output_angle = 0


# angle_out = []
m_max = 1000000000

