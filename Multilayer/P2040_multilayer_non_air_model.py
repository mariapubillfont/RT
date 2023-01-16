# Test the P2040 multilayer model for non-air material in last layer

import numpy as np
import matplotlib.pyplot as plt
import P2040_model as model

# General parameters
frequency = 13e9
lambd = 3e8/frequency
k_0 = 2*np.pi/lambd
d_core = 0.1 # Core thickness
er = 2.5 # Core permittivity



## First matching interface
incidenceAngleRadians1 = np.linspace(0,np.pi/2,91)
layerThickness = [0, lambd/(4*np.sqrt(np.sqrt(er))), 0] # Matching layer
complexRelativePermittivity = [1, np.sqrt(er), er] # From air to dielectric
rTE1 = np.ones([np.size(incidenceAngleRadians1),1])
tTE1 = np.ones([np.size(incidenceAngleRadians1),1])
rTM1 = np.ones([np.size(incidenceAngleRadians1),1])
tTM1 = np.ones([np.size(incidenceAngleRadians1),1])
for k in range(0,np.size(incidenceAngleRadians1)):
    rTE1[k],rTM1[k],tTE1[k],tTM1[k] = model.P2040multilayer(k_0,incidenceAngleRadians1[k],layerThickness,complexRelativePermittivity)




## Core propagation
# This is different for a plane wave and a ray?
#
# For a plane wave, the reference for the propagation phase is normal to
# the parallel surfaces of the core slab, and in the limit of grazing, the
# phase difference between upper and lower surface is zero for a plane
# wave. Hence, the np.expression for phase delay then scales with np.cos(theta).
#
# For a ray, the travel distance is related to the path length travelled in
# the core. Hence, the phase delay scales with 1/np.cos(theta), since the ray
# is travelling longer than the thickness of the core for non-normal
# (theta>0) incidence.
incidenceAngleRadians2 = np.asin(np.sin(incidenceAngleRadians1)/np.sqrt(er))
tTE2 = np.exp(1j * k_0 * d_core / np.cos(incidenceAngleRadians2))
tTM2 = np.exp(1j * k_0 * d_core / np.cos(incidenceAngleRadians2))



## Second matching interface
incidenceAngleRadians3 = incidenceAngleRadians1
layerThickness = [0, lambd/(4*np.sqrt(np.sqrt(er))), 0] # Matching layer
complexRelativePermittivity = [er, np.sqrt(er), 1] # From dielectric to air
rTE3 = np.ones(np.size(incidenceAngleRadians3),1)
tTE3 = np.ones(np.size(incidenceAngleRadians3),1)
rTM3 = np.ones(np.size(incidenceAngleRadians3),1)
tTM3 = np.ones(np.size(incidenceAngleRadians3),1)
for k in range(1,np.size(incidenceAngleRadians1)):
    [rTE3[k],rTM3[k],tTE3[k],tTM3[k]] = model.P2040multilayer(k_0,incidenceAngleRadians3[k],layerThickness,complexRelativePermittivity)



## Power conservation tests
# Compare (34)-(36) in https://www.itu.int/rec/R-REC-P.2040-2-202109-I
#
# I don't fully understand this result. Why is the rescaling the same for
# both entering and leaving the core dielectric. I thought this was related
# to the E-field scaling with the dielectric constant... Hmm, D-field?



# First matching layer rescaling for power conservation test
projScale = np.cos(incidenceAngleRadians2)/np.cos(incidenceAngleRadians1)
rescaleDiel = np.sqrt(er/1) # From air to dielectric
plt.figure(1)
plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.np.log10(abs(rTE1)),'k-')
plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(tTE1)))
plt.plot(np.rad2deg(incidenceAngleRadians1), 10*np.log10(abs(rTE1)*abs(rTE1) + 1/rescaleDiel*projScale*abs(tTE1)*abs(tTE1)))
plt.grid()
plt.ylim([-15, 5.1])
plt.show()



plt.figure(2)
plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(rTM1)),'k-')
plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(tTM1)))
plt.plot(np.rad2deg(incidenceAngleRadians1), 10*np.log10(abs(rTM1)*abs(rTM1) +  rescaleDiel*projScale*abs(tTM1)*abs(tTM1)))
plt.grid()
plt.ylim([-15, 5.1])
plt.show()


# Second matching layer rescaling for power conservation test
projScale = np.cos(incidenceAngleRadians3)/np.cos(incidenceAngleRadians2)
rescaleDiel = np.sqrt(er/1) # From dielectric to air but this is the same?!
plt.figure(3)
plt.plot(np.rad2deg(incidenceAngleRadians2), 20*np.log10(abs(rTE3)),'k-')
plt.plot(np.rad2deg(incidenceAngleRadians2), 20*np.log10(abs(tTE3)))
plt.plot(np.rad2deg(incidenceAngleRadians2), 10*np.log10(abs(rTE3)*abs(rTE3) + 1/rescaleDiel*projScale*abs(tTE3)*abs(tTE3)))
plt.grid()
plt.ylim([-15, 5.1])
plt.show()


plt.figure(4)
plt.plot(np.rad2deg(incidenceAngleRadians2), 20*np.log10(abs(rTM3)),'k-')
plt.plot(np.rad2deg(incidenceAngleRadians2), 20*np.log10(abs(tTM3)))
plt.plot(np.rad2deg(incidenceAngleRadians2), 10*np.log10(abs(rTM3)*abs(rTM3) +  rescaleDiel*projScale*abs(tTM3)*abs(tTM3)))
plt.grid()
plt.ylim([-15, 5.1])
plt.show()


## Transmitted power, air-dielectric-air
# I do not understand how er comes in here. The different factors for TE
# and TM are simply added to get transmitted power to 1 at normal
# incidence. Why would they differ for TE and TM??! How does this relate to
# the continuity of the tangential E-field and normal D-field?
plt.figure(5)
plt.plot(np.rad2deg(incidenceAngleRadians1), 10*np.log10(abs(tTE1 * tTE2 * tTE3)^2),'k-')
plt.plot(np.rad2deg(incidenceAngleRadians1), 10*np.log10(abs(tTM1 * tTM2 * tTM3)^2),'r--')
plt.grid()
plt.ylim([-15, 5.1])
plt.show()
