
import numpy as np
import matplotlib.pyplot as plt
f = 13
c0 = 0.299792458 
wv = c0/f
k_0 = 2*np.pi/(c0/f)
a_in = np.deg2rad(0)
#thickness = np.arange(0,0.1,0.00001)
thickness = 0.1
#n2 = 5.31
c = 0.0326
d = 0.8095
cond = c*f**d #S/m
er = 25
#permittivity = n2 - 1j*(17.98*cond/f)
#tan_delta = 0.00066
tan_delta = 0
permittivity = er*(1-1j*tan_delta)

#thickness = 0.02
#freq = np.arange(1, 20, 0.01) #freq in GHz
#wv_range = c0/freq




def getReflectionCoefficients(wv, thickness, polaritzation, permittivity, incidentAngle):
    eta = np.sqrt(permittivity-np.sin(incidentAngle)**2)
    if polaritzation == 'TE': R_aux = (np.cos(incidentAngle)-eta)/(np.cos(incidentAngle)+eta)
    else: R_aux = (permittivity*np.cos(a_in)-eta)/(permittivity*np.cos(incidentAngle)+eta)
    q = eta*2*np.pi*thickness/wv
    T = ((1-R_aux**2)*np.exp(-1j*q))/(1-R_aux**2*np.exp(-2j*q))
    R = R_aux*(1-np.exp(-2j*q))/(1-R_aux**2*np.exp(-2j*q))
    return T, R

mu_0 = np.pi*4e-7
mu_r1 = 1
mu_r2 = 1
ep_0 = 8.854e-12
ep_r1 = 1
ep_r2 = 25
z1 = np.sqrt((mu_0)/(ep_0))
z2 = np.sqrt((mu_r2)/(ep_r2))
reflection = (z2-z1)/(z2+z1)


T, R = getReflectionCoefficients(wv, thickness, 'TE', permittivity, a_in)
plt.plot(thickness, 20*np.log10(abs(R)))
plt.plot(thickness, 20*np.log10(abs(T)))
plt.xlabel('Slab thickness (m)')
print(20*np.log10(abs(R)), 20*np.log10(abs(T)))
# plt.plot(c0/wv_range, 20*np.log10(abs(R)))
# plt.plot(c0/wv_range, 20*np.log10(abs(T)))
#plt.xlabel('Frequency (GHz)')
plt.ylabel('T and R coefficents amplitude (dB)')

plt.title('Transmission and reflection coefficents')
plt.legend(['Reflection', 'Tramission'])

plt.grid()
plt.show()
