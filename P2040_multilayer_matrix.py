
# Test the P2040 multilayer model for non-air material in last layer

import numpy as np
import matplotlib.pyplot as plt
#import P2040_model_matrix as modmat
import input as I
#import model_multilayer as multilayer


# General parameters
# frequency = 3e8
# lambd = 3e8/frequency
# k_0 = 2*np.pi/lambd
# d_core = lambd # Core thickness
# tan_delta = 0.0
# er = 2.5*(1 - 1j*tan_delta) # Core permittivity

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
# Transfer matrix for multilayer slab consisting of N finite thickness
# layers, for N+2 permittivities (for the initial region, the N layers, and
# the final region).
#
# The transfer matrix is refered to the interface between layer 0 (the
# incident region) and layer 1, and the interface between layer N and layer
# N+1 (the final region). The thickness values of layers 0 and N+1 are
# hence assumed to be zero (and not used).
#



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



def getReflectionCoefficients_ML(incidentAngle, layerThickness_in, er, frequency):
    lambd = 299792458/frequency
    k_0 = 2*np.pi/lambd
    d_core = layerThickness_in[1]
    ################################# First matching interface #####################################################################
    #incidenceAngleRadians1 = np.transpose(np.linspace(0,np.pi/2,91))
    incidenceAngleRadians1 = incidentAngle[0]
    layerThickness = [0, layerThickness_in[0], 0] # Matching layer
    complexRelativePermittivity = [1, np.sqrt(er), er] # From air to dielectric

    # Transfer matrix model
    rTE1 = np.ones([np.size(incidenceAngleRadians1),1],dtype=np.complex_)
    tTE1 = np.ones([np.size(incidenceAngleRadians1),1], dtype=np.complex_)
    rTM1 = np.ones([np.size(incidenceAngleRadians1),1],dtype=np.complex_)
    tTM1 = np.ones([np.size(incidenceAngleRadians1),1], dtype=np.complex_)
    A = multiLayerTransferMatrix(incidenceAngleRadians1, layerThickness, complexRelativePermittivity, frequency, 'te')
    rTE1 = A[1][0]/A[0][0]
    tTE1 = 1/A[0][0]
    A = multiLayerTransferMatrix(incidenceAngleRadians1, layerThickness, complexRelativePermittivity, frequency, 'tm')
    rTM1 = A[1,0]/A[0,0]
    tTM1 = 1/A[0,0]

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
    tTE2 = np.exp(-1j * k_0*np.sqrt(er) * d_core)
    tTM2 = np.exp(1j * k_0/np.sqrt(er) * d_core)
    #############################################################################################################################





    ################################### Second matching interface#################################################################
    layerThickness = [0,  layerThickness_in[2], 0] # Matching layer
    complexRelativePermittivity = [er, np.sqrt(er), 1] # From dielectric to air
    
    incidenceAngleRadians2 = incidentAngle[2]
    # Transfer matrix model
    rTE3 = np.ones([np.size(incidenceAngleRadians2),1], dtype=np.complex_)
    tTE3 = np.ones([np.size(incidenceAngleRadians2),1], dtype=np.complex_)
    rTM3 = np.ones([np.size(incidenceAngleRadians2),1], dtype=np.complex_)
    tTM3 = np.ones([np.size(incidenceAngleRadians2),1],dtype=np.complex_)
    A = multiLayerTransferMatrix(incidenceAngleRadians2, layerThickness, complexRelativePermittivity, frequency, 'te')
    rTE3 = A[1,0]/A[0,0]
    tTE3 = 1/A[0,0]
    A = multiLayerTransferMatrix(incidenceAngleRadians2, layerThickness, complexRelativePermittivity, frequency, 'tm')
    rTM3 = A[1,0]/A[0,0]
    tTM3 = 1/A[0,0]
    #############################################################################################################################


    return (tTE1* tTE2 * tTE3)