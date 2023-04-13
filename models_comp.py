 
# Test the P2040 multilayer model for non-air material in last layer

import numpy as np
import matplotlib.pyplot as plt
#import P2040_model_matrix as modmat
import input as I
#import model_multilayer as multilayer


def fresnelCoefficients_TE(theta1, theta2, n1, n2):
    tmp1 = n1 * np.cos(theta1)
    tmp2 = n2 * np.cos(theta2)
    r_te = (tmp1 - tmp2) / (tmp1 + tmp2)
    t_te = 2 * tmp1 / (tmp1 + tmp2)
    return [r_te, t_te] 

def fresnelCoefficients_TM(theta1, theta2, n1, n2):
    tmp1 = n1 * np.cos(theta2)
    tmp2 = n2 * np.cos(theta1)
    r_te = (tmp1 - tmp2) / (tmp1 + tmp2)
    t_te = 2 * n1 * np.cos(theta1) / (tmp1 + tmp2)
    return [r_te, t_te]


def multiLayerTransferMatrix(theta_in, t, er, frequency, pol):
    # Wavelength
    lambd = 299792458/frequency
    k0 = 2*np.pi/lambd

    # Transfer matrix
    A = np.identity(2)
    for n in range (1,np.size(er)):
        # Angle in outgoing layer (5.4)
        theta_out = np.arcsin(np.sqrt(er[n-1]/er[n])*np.sin(theta_in))

    # Propagation constant normal to boundary ok for plane wave for ray?!
        k_n = k0 * np.sqrt(er[n]) * np.cos(theta_out)

    # Fresnel reflection and transmission coefficients
        if pol == 'te':
            [R, T] = fresnelCoefficients_TE(theta_in, theta_out, np.sqrt(er[n-1]), np.sqrt(er[n]))
        else:
            [R, T] = fresnelCoefficients_TM(theta_in, theta_out, np.sqrt(er[n-1]), np.sqrt(er[n]))
        
    # Transfer matrix
        e = 1j * k_n * t[n]
        A = 1/T  * np.matmul(A, [ [np.exp(e), R*np.exp(-e)],  [R*np.exp(e),  np.exp(-e)]])

    # Update angle
        theta_in = theta_out

    return A



def getReflectionCoefficients_TMM(incidentAngle, layerThickness_in, complexRelativePermittivity, frequency):
    incidenceAngleRadians1 = np.deg2rad(incidentAngle)
    layerThickness = layerThickness_in
  
    # Transfer matrix model
    rTE1 = np.ones([np.size(incidenceAngleRadians1),1],dtype=np.complex_)
    tTE1 = np.ones([np.size(incidenceAngleRadians1),1], dtype=np.complex_)
    rTM1 = np.ones([np.size(incidenceAngleRadians1),1],dtype=np.complex_)
    tTM1 = np.ones([np.size(incidenceAngleRadians1),1], dtype=np.complex_)
    for i in range(0, len(incidentAngle)):
        A = multiLayerTransferMatrix(incidenceAngleRadians1[i], layerThickness, complexRelativePermittivity, frequency, 'te')
        rTE1[i] = A[1][0]/A[0][0]
        tTE1[i] = 1/A[0][0]
        #tTE1[i] = 1/A[0][0]*np.exp(-1j*k0*(I.thickness_ML1*2+lambd/np.sqrt(er))*np.tan(incidenceAngleRadians1[i])*np.sin(incidenceAngleRadians1[i]))
    return tTE1
 


def getReflectionCoefficients_ITU(k_0, layerThickness, polaritzation, complexPermittivity, incidentAngle):
   # layerThickness = [0,layerThicknessi, 0 ]
    nLayers = len(layerThickness)                                                           #number of layers
    eta_n = complexPermittivity                                                             #assuming permeability mu_0
    k_n = k_0*np.sqrt(eta_n)                                                                #wave length in material
    incidenceAngle_n = np.arcsin(np.sin(incidentAngle)/np.sqrt(eta_n))                      #
    N = nLayers

    # % Initialize
    A =  np.ones(N, dtype=np.complex_)
    B = np.zeros(N, dtype=np.complex_)
    F =  np.ones(N, dtype=np.complex_)
    G = np.zeros(N, dtype=np.complex_)
    
##Calculate backwards from last layer to first
    for i in range(N-2,-1,-1):
        W = np.cos(incidenceAngle_n[i+1])/np.cos(incidenceAngle_n[i])*np.sqrt(eta_n[i]/eta_n[i+1])
        Y = np.cos(incidenceAngle_n[i+1])/np.cos(incidenceAngle_n[i])*np.sqrt(eta_n[i+1]/eta_n[i])
        A[i] = 0.5*np.exp(1j*k_n[i]*layerThickness[i]*np.cos(incidenceAngle_n[i]))*(A[i+1]*(1+Y)+B[i+1]*(1-Y))
        B[i] = 0.5*np.exp(-1j*k_n[i]*layerThickness[i]*np.cos(incidenceAngle_n[i]))*(A[i+1]*(1-Y)+B[i+1]*(1+Y))
        F[i] = 0.5*np.exp( 1j*k_n[i]*layerThickness[i]*np.cos(incidenceAngle_n[i]))*(F[i+1]*(1+W)+G[i+1]*(1-W))
        G[i] = 0.5*np.exp(-1j*k_n[i]*layerThickness[i]*np.cos(incidenceAngle_n[i]))*(F[i+1]*(1-W)+G[i+1]*(1+W))

    if polaritzation == 'TE':
        r = B[0]/A[0] #reflection and transmission for TE polaritzation
        t = 1/A[0]
    else:
        r = G[0]/F[0] #reflection and transmission for TM polartization
        t = 1/F[0]
    return t





incidentAngle = np.linspace(0,80,80)
#incidentAngle = 60
#thickness = [0, I.thickness_ML1, 0.1, I.thickness_ML1, 0]
f = 13e9
lambd = 299792458/f
k0 = 2*np.pi/lambd
er=2.5

#thickness = [0, I.thickness_ML1, lambd/np.sqrt(er), I.thickness_ML1,0]
thickness = [0, I.thickness_ML1, lambd/np.sqrt(er)]
cond = 0.0326*1**0.8095
complex_er = er - 1j*17.98*cond/1
complexRelativePermittivity = [1, np.sqrt(er), er, np.sqrt(er), 1] # From air to dielectric
complexRelativePermittivity = [1, np.sqrt(er), er]




#thickness = [0, lambd/(2*np.sqrt(er)), 0]
#thickness = np.linspace(0.0001,1,1000)
T_ITU = []
T_TMM = getReflectionCoefficients_TMM(incidentAngle, thickness, complexRelativePermittivity, f)
for i in range(0, len(incidentAngle)):
    T_ITU = np.append(T_ITU, getReflectionCoefficients_ITU(k0, thickness, 'TE', complexRelativePermittivity, np.deg2rad(incidentAngle[i])))

plt.figure()
plt.plot(incidentAngle, np.abs(T_TMM))
plt.plot(incidentAngle, np.abs(T_ITU))
#plt.ylim([-30, 0])
plt.legend(("abs(TMM)", "abs(ITU)"))
plt.grid()
plt.show()


# plt.figure()
# # plt.plot(incidentAngle, np.abs(T_TMM))
# plt.plot(thickness, np.angle(T_ITU))
# plt.legend(("ph(TMM)", "abs(ITU)"))
# plt.grid()
# plt.show()
plt.figure()
plt.plot(incidentAngle, np.angle(T_TMM))
plt.plot(incidentAngle, np.angle(T_ITU))
plt.legend(('ph(TMM)', 'ph(ITU)'))
plt.grid()
plt.show()