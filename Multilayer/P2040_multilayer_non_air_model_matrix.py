
# Test the P2040 multilayer model for non-air material in last layer

import numpy as np
import matplotlib.pyplot as plt
import reflections_matrix_model as modmat
#import model_multilayer as multilayer


# General parameters
frequency = 13e9
lambd = 3e8/frequency
k_0 = 2*np.pi/lambd
d_core = 0.1 # Core thickness
tan_delta = 0
er = 25*(1 - 1j*tan_delta) # Core permittivity




def getReflectionCoefficients_multiLayer(k_0, layerThickness, polaritzation, complexPermittivity, incidentAngle):
    nLayers = len(layerThickness)
    eta_n = complexPermittivity
    k_n = k_0*np.sqrt(eta_n)
    rTE = np.zeros(len(incidentAngle), dtype=np.complex_)
    tTE = np.zeros(len(incidentAngle), dtype=np.complex_)
    rTM = np.zeros(len(incidentAngle), dtype=np.complex_)
    tTM = np.zeros(len(incidentAngle), dtype=np.complex_)
    for j in range(0, len(incidentAngle)):

        incidenceAngle_n = np.arcsin(np.sin(incidentAngle[j])/np.sqrt(eta_n))
        N = nLayers

        # % Initialize
        A =  np.ones(N, dtype=np.complex_)##*np.sqrt(complexPermittivity[N-1])
        B = np.zeros(N, dtype=np.complex_)
        F =  np.ones(N, dtype=np.complex_)##/np.sqrt(complexPermittivity[N-1])
        G = np.zeros(N, dtype=np.complex_)
        
    ##Calculate backwards from last layer to first
        for i in range(N-2,-1,-1):
            W = np.cos(incidenceAngle_n[i+1])/np.cos(incidenceAngle_n[i])*np.sqrt(eta_n[i]/eta_n[i+1])
            Y = np.cos(incidenceAngle_n[i+1])/np.cos(incidenceAngle_n[i])*np.sqrt(eta_n[i+1]/eta_n[i])
            A[i] = 0.5*np.exp(1j*k_n[i]*layerThickness[i]*np.cos(incidenceAngle_n[i]))*(A[i+1]*(1+Y)+B[i+1]*(1-Y))
            B[i] = 0.5*np.exp(-1j*k_n[i]*layerThickness[i]*np.cos(incidenceAngle_n[i]))*(A[i+1]*(1-Y)+B[i+1]*(1+Y))
            F[i] = 0.5*np.exp( 1j*k_n[i]*layerThickness[i]*np.cos(incidenceAngle_n[i]))*(F[i+1]*(1+W)+G[i+1]*(1-W))
            G[i] = 0.5*np.exp(-1j*k_n[i]*layerThickness[i]*np.cos(incidenceAngle_n[i]))*(F[i+1]*(1-W)+G[i+1]*(1+W))

        rTE[j] = B[0]/A[0]
        rTM[j] = G[0]/F[0]
        tTE[j] = 1/A[0]
        tTM[j] = 1/F[0]
    if polaritzation == 'TE' : 
        return tTE, rTE
    else:
        return tTM, rTM    




############################### multilayer model #####################################################################
incidenceAngleRadians = np.transpose(np.linspace(0,np.pi/2,91))
layerThickness = [0, lambd/(4*np.sqrt(np.sqrt(abs(er)))), d_core, lambd/(4*np.sqrt(np.sqrt(abs(er)))), 0] # Matching layer
complexRelativePermittivity = [1, np.sqrt(er), er, np.sqrt(er), 1] # From air to dielectric

# layerThickness = [0, d_core, 0]
# complexRelativePermittivity = [1, er, 1]


tTM_multilayer = np.zeros(len(incidenceAngleRadians), dtype=np.complex_)
rTM_multilayer = np.zeros(len(incidenceAngleRadians), dtype=np.complex_)
tTE_multilayer = np.zeros(len(incidenceAngleRadians), dtype=np.complex_)
rTE_multilayer = np.zeros(len(incidenceAngleRadians), dtype=np.complex_)
tTE_multilayer, rTE_multilayer = getReflectionCoefficients_multiLayer(k_0, layerThickness, 'TE', complexRelativePermittivity, incidenceAngleRadians)
# plt.figure()
# plt.ylim([-20, 1.1])
# plt.plot(np.rad2deg(incidenceAngleRadians),  20*np.log10(abs(rTE_multilayer)))
# plt.plot(np.rad2deg(incidenceAngleRadians),  20*np.log10(abs(tTE_multilayer)))
# plt.plot(np.rad2deg(incidenceAngleRadians), 10*np.log10(abs(tTE_multilayer)**2 + abs(rTE_multilayer)**2))
# plt.legend(['Reflection', 'Transmission'])
# plt.grid()
# plt.show()

tTM_multilayer, rTM_multilayer = getReflectionCoefficients_multiLayer(k_0, layerThickness, 'TM', complexRelativePermittivity, incidenceAngleRadians)

# #########################################################################################################################



################################# First matching interface #####################################################################
#incidenceAngleRadians1 = np.transpose(np.linspace(0,np.pi/2,91))
incidenceAngleRadians1 =  np.transpose(np.linspace(0,np.pi/2,91))
layerThickness = [0, lambd/(4*np.sqrt(np.sqrt(abs(er)))), 0] # Matching layer
complexRelativePermittivity = [1, np.sqrt(er), er] # From air to dielectric

# layerThickness = [0, 0 ]
# complexRelativePermittivity = [1, er]

# Transfer matrix model
rTE1 = np.ones([np.size(incidenceAngleRadians1),1],dtype=np.complex_)
tTE1 = np.ones([np.size(incidenceAngleRadians1),1], dtype=np.complex_)
rTM1 = np.ones([np.size(incidenceAngleRadians1),1],dtype=np.complex_)
tTM1 = np.ones([np.size(incidenceAngleRadians1),1], dtype=np.complex_)
for k in range(0, np.size(incidenceAngleRadians1)):
    A = modmat.multiLayerTransferMatrix(incidenceAngleRadians1[k], layerThickness, complexRelativePermittivity, frequency, 'te')
    rTE1[k] = A[1][0]/A[0][0]
    tTE1[k] = 1/A[0][0]
    A = modmat.multiLayerTransferMatrix(incidenceAngleRadians1[k], layerThickness, complexRelativePermittivity, frequency, 'tm')
    rTM1[k] = A[1,0]/A[0,0]
    tTM1[k] = 1/A[0,0]

#############################################################################################################################



####################################### Core propagation #####################################################################
# Is this different for a plane wave and a ray?!
#
# For a plane wave, the reference for the propagation phase is normal to
# the parallel surfaces of the core slab, and in the limit of grazing, the
# phase difference between upper and lower surface is zero for a plane
# wave. Hence, the np.expression for phase delay then scales with cos(theta).
#
# For a ray, the travel distance is related to the path length travelled in
# the core. Hence, the phase delay scales with 1/cos(theta), since the ray
# is travelling longer than the thickness of the core for oblique
# (theta>0) incidence.
incidenceAngleRadians2 = np.arcsin(np.sin(incidenceAngleRadians1)/np.sqrt(er))
#tTE2 = np.exp(-1j * k_0*np.sqrt(er) * d_core * np.cos(incidenceAngleRadians2)) 
tTE2 = np.exp(-1j * k_0*np.sqrt(er) * d_core) 
#treure el cosine i posar d2 directment, mirar la phase aver si surten els mateixos valors
tTM2 = np.exp(-1j * k_0*np.sqrt(er) * d_core)
#############################################################################################################################





################################### Second matching interface#################################################################
# layerThickness = [0, 0] # Matching layer
# complexRelativePermittivity = [er, 1] # From dielectric to air

layerThickness = [0, lambd/(4*np.sqrt(np.sqrt(abs(er)))), 0] # Matching layer
complexRelativePermittivity = [er, np.sqrt(er), 1] # From air to dielectric

# Transfer matrix model
rTE3 = np.ones([np.size(incidenceAngleRadians2),1], dtype=np.complex_)
tTE3 = np.ones([np.size(incidenceAngleRadians2),1], dtype=np.complex_)
rTM3 = np.ones([np.size(incidenceAngleRadians2),1], dtype=np.complex_)
tTM3 = np.ones([np.size(incidenceAngleRadians2),1],dtype=np.complex_)
for k in range(0, np.size(incidenceAngleRadians1)):
    A = modmat.multiLayerTransferMatrix(incidenceAngleRadians2[k], layerThickness, complexRelativePermittivity, frequency, 'te')
    rTE3[k] = A[1,0]/A[0,0]
    tTE3[k] = 1/A[0,0]
    A = modmat.multiLayerTransferMatrix(incidenceAngleRadians2[k], layerThickness, complexRelativePermittivity, frequency, 'tm')
    rTM3[k] = A[1,0]/A[0,0]
    tTM3[k] = 1/A[0,0]
#############################################################################################################################





# First matching layer rescaling for power conservation test
incidenceAngleRadians3 = incidenceAngleRadians1
projScale = np.cos(incidenceAngleRadians2)/np.cos(incidenceAngleRadians3)
rescaleDiel = np.sqrt(er/1) # From air to dielectric
# plt.figure(1)
# plt.clf()
# plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(rTE1)),'k-', linewidth= 1)
# plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(tTE1)),'r-', linewidth= 1)
# plt.plot(np.rad2deg(incidenceAngleRadians1), 10*np.log10((abs(rTE1)*abs(rTE1))[:,0] + np.multiply((abs(rescaleDiel*projScale)),(abs(tTE1)*abs(tTE1))[:,0])), 'b-', linewidth=  1)
# plt.ylim([-15.1, 10.1])
# plt.title('First matching layer: 1 -> np.sqrt(er) -> er')
# plt.legend(['R (TE)', 'T (TE)', 'P (TE)'], loc= 'lower left')
# plt.grid()
# plt.show()

# plt.figure(2)
# plt.clf()
# plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(rTM1)),'k-', linewidth=  1)
# plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(tTM1)),'r-', linewidth=  1)
# plt.plot(np.rad2deg(incidenceAngleRadians1), 10*np.log10((abs(rTM1)*abs(rTM1))[:,0] + np.multiply((abs(rescaleDiel*projScale)),(abs(tTM1)*abs(tTM1))[:,0])), 'b-', linewidth=  1)
# plt.ylim([-15.1, 10.1])
# plt.title('First matching layer: 1 -> np.sqrt(er) -> er')
# plt.legend(['R (TM)', 'T (TM)', 'P (TM)'], loc= 'lower left')
# plt.grid()
# #plt.show()

# # Second matching layer rescaling for power conservation test
# projScale = np.cos(incidenceAngleRadians3)/np.cos(incidenceAngleRadians2)
# rescaleDiel = np.sqrt(1/er) # From dielectric to air but this is the same?!
# plt.figure(3)
# plt.clf()
# plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(rTE3)),'k-', linewidth=  1)
# plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(tTE3)),'r-', linewidth=  1)
# plt.plot(np.rad2deg(incidenceAngleRadians1), 10*np.log10((abs(rTE3)*abs(rTE3))[:,0] + np.multiply((abs(rescaleDiel*projScale)),(abs(tTE3)*abs(tTE3))[:,0])), 'b-', linewidth=  1)
# plt.ylim([-15.1, 10.1])
# plt.title('Second matching layer: er -> np.sqrt(er) -> 1')
# plt.legend(['R (TE)', 'T (TE)', 'P (TE)'], loc= 'lower left')
# plt.grid()
# #plt.show()

# plt.figure(4)
# plt.clf()
# plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(rTM3)),'k-', linewidth=  1)
# plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(tTM3)),'r-', linewidth=  1)
# plt.plot(np.rad2deg(incidenceAngleRadians1), 10*np.log10((abs(rTM3)*abs(rTM3))[:,0] + np.multiply((abs(rescaleDiel*projScale)),(abs(tTM3)*abs(tTM3))[:,0])), 'b-', linewidth=  1)
# plt.ylim([-15.1, 10.1])
# plt.title('Second matching layer: er -> np.sqrt(er) -> 1')
# plt.legend(['R (TM)', 'T (TM)', 'P (TM)'], loc= 'lower left')
# plt.grid()
#plt.show()


## Transmitted power, air-dielectric-air
# I do not understand how er comes in here. The different factors for TE
# and TM are simply added to get transmitted power to 1 at normal
# incidence. Why would they differ for TE and TM??! How does this relate to
# the continuity of the tangential E-field and normal D-field?
fig = plt.figure(5)
plt.clf()
plt.rcParams["font.family"] = "Times New Roman" 
ax1 = fig.add_subplot(111)

ax1.xaxis.label.set_fontsize(14)
ax1.yaxis.label.set_fontsize(14)
plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(tTE_multilayer)), color = '#0072BD', linewidth = 2 )
plt.plot(np.rad2deg(incidenceAngleRadians1), 10*np.log10(abs(tTE1[:,0] * tTE2 * tTE3[:,0])**2), color = '#77AC30', linewidth = 2)
#plt.plot(np.rad2deg(incidenceAngleRadians1), -3*np.ones(np.size(incidenceAngleRadians1)), 'g-', linewidth = 1 )
plt.plot(np.rad2deg(incidenceAngleRadians1), 20*np.log10(abs(tTM_multilayer)), ':',color = '#0072BD', linewidth = 2 )
plt.plot(np.rad2deg(incidenceAngleRadians1), 10*np.log10(abs(tTM1[:,0] * tTM2 * tTM3[:,0])**2),':',color = '#77AC30', linewidth = 2)

plt.ylim([-15.1, 0])
plt.title('ML in Slab of thickness = 100 mm and $\epsilon_r$ = 25')
#plt.legend(['T (TE)', 'T (TE) P20140 model','T (TM)', 'T (TM) multilayer'], loc= 'lower left')
plt.legend(['T (TE) P20140 ITU model', 'T (TE) double interface model', 'T (TM) P20140 ITU model', 'T (TM) double interface model' ], loc= 'lower left')
plt.ylabel('dB')
plt.xlabel('Incident angle (deg)')

plt.grid()
plt.show()
