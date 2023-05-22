 
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

def multiLayerTransferMatrixMod(theta_in, t, er, frequency, pol):
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
        k_n = k0 * np.sqrt(er[n]) * np.cos(theta_out) #SHOULD BE A *COS

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
        k_n = k0 * np.sqrt(er[n]) * np.cos(theta_out) #SHOULD BE A *COS

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


################################# AGGREGATE #####################################################################
def getReflectionCoefficients_agg(incidentAngle,layerThickness_in,complexRelativePermittivity, frequency, delta):
    lambd = 299792458/frequency
    k_0 = 2*np.pi/lambd
    if len(layerThickness_in) == 1:
        layerThickness = np.pad(layerThickness_in, (1,1), 'constant', constant_values=(0,0))
        incidentAngle = np.pad(incidentAngle, (1,1), 'edge')
        distance = delta[0]
    else:
        layerThickness = [0, layerThickness_in[0], layerThickness_in[1], layerThickness_in[2], 0]    
        distance = delta[2] 
    incidenceAngleRadians1 = incidentAngle[0]
    #layerThickness = [0, layerThickness_in[0], layerThickness_in[1], layerThickness_in[2], 0]
    #complexRelativePermittivity = [1, np.sqrt(er), er, np.sqrt(er), 1] # From air to dielectric
    #theta_out =np.arcsin(np.sqrt(np.sqrt(er))*np.sin(incidentAngle[3]))
    #ni/nout*sin(oi)
    # Transfer matrix model
    rTE1 = np.ones([np.size(incidenceAngleRadians1),1],dtype=np.complex_)
    tTE1 = np.ones([np.size(incidenceAngleRadians1),1], dtype=np.complex_)
    rTM1 = np.ones([np.size(incidenceAngleRadians1),1],dtype=np.complex_)
    tTM1 = np.ones([np.size(incidenceAngleRadians1),1], dtype=np.complex_)
    A = multiLayerTransferMatrixMod(incidenceAngleRadians1, layerThickness, complexRelativePermittivity, frequency, 'te')
    rTE1 = A[1][0]/A[0][0]
    tTE1 = 1/A[0][0]*np.exp(-1j*k_0*distance*np.sin(abs(incidenceAngleRadians1)))
    A = multiLayerTransferMatrix(incidenceAngleRadians1, layerThickness, complexRelativePermittivity, frequency, 'tm')
    rTM1 = A[1,0]/A[0,0]
    tTM1 = 1/A[0,0]
    return tTE1




################################# CASCADE #####################################################################
def getReflectionCoefficients_cascade(incidentAngle,layerThickness_in,er, frequency):
    phase_compensation = 0
    lambd = 299792458/frequency
    k_0 = 2*np.pi/lambd
    if len(layerThickness_in) == 1:
        d_core = layerThickness_in[0]
       # layerThickness_in = np.pad(layerThickness_in, (1,1), 'constant', constant_values=(0,0))
        incidentAngle = np.pad(incidentAngle, (1,1), 'edge')
    else:
        d_core = layerThickness_in[1]
    ################################# First matching interface #####################################################################
    incidenceAngleRadians1 = incidentAngle[0]

    if len(layerThickness_in) == 1:
        layerThickness = [0, 0]
        complexRelativePermittivity = [1, er]
        phase_compensation = 1
    else:
        layerThickness = [0, layerThickness_in[0], 0]
        complexRelativePermittivity = [1, np.sqrt(er), (er)]
        phase_compensation = 1
    
    theta_out = np.arcsin(np.sin(incidentAngle[0])/np.sqrt(np.sqrt(er))) #refracted angle inisde the matching layer
    #Z_diec= np.sqrt(2.5-np.sin(incidenceAngleRadians1)**2)/(2.5*np.cos(incidenceAngleRadians1))
    # Transfer matrix model
    rTE1 = np.ones([np.size(incidenceAngleRadians1),1],dtype=np.complex_)
    tTE1 = np.ones([np.size(incidenceAngleRadians1),1], dtype=np.complex_)
    rTM1 = np.ones([np.size(incidenceAngleRadians1),1],dtype=np.complex_)
    tTM1 = np.ones([np.size(incidenceAngleRadians1),1], dtype=np.complex_)
    A = multiLayerTransferMatrix(incidenceAngleRadians1, layerThickness, complexRelativePermittivity, frequency, 'te')
    rTE1 = A[1][0]/A[0][0]
    if phase_compensation == 1:
        tTE1 = 1/A[0][0]*np.exp(-1j*k_0*(I.thickness_ML1)*np.tan(theta_out)*np.sin(theta_out))
        #*np.sqrt(Z_diec)
    else:
        tTE1 = 1/A[0][0]   
    #/(Z_diec)
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
    #incidenceAngleRadians2 = np.arcsin(np.sin(incidenceAngleRadians1)/np.sqrt(er))
    tTE2 = np.exp(-1j * k_0*np.sqrt(er) * d_core)
    tTM2 = np.exp(1j * k_0*np.sqrt(er) * d_core)
    #############################################################################################################################

    ################################### Second matching interface#################################################################
    #layerThickness = [0,layerThickness_in[2], 0] # Matching layer
    if len(layerThickness_in) == 1:
        layerThickness = [d_core, 0]
        complexRelativePermittivity = [er, 1]
        phase_compensation = 0
    else:
        complexRelativePermittivity = [er,np.sqrt(er), 1] # From dielectric to air
        layerThickness = [0, layerThickness_in[2], 0]
        phase_compensation = 1


    
    incidenceAngleRadians2 = incidentAngle[2]
    theta_out2 = np.arcsin(np.sqrt(np.sqrt(er))*np.sin(incidentAngle[2]))
    Z_diec2 = np.sqrt(1-np.sin(incidenceAngleRadians2)**2)/(1*np.cos(incidenceAngleRadians2))

    # Transfer matrix model
    rTE3 = np.ones([np.size(incidenceAngleRadians2),1], dtype=np.complex_)
    tTE3 = np.ones([np.size(incidenceAngleRadians2),1], dtype=np.complex_)
    rTM3 = np.ones([np.size(incidenceAngleRadians2),1], dtype=np.complex_)
    tTM3 = np.ones([np.size(incidenceAngleRadians2),1],dtype=np.complex_)
    A = multiLayerTransferMatrix(incidenceAngleRadians2, layerThickness, complexRelativePermittivity, frequency, 'te')
    rTE3 = A[1,0]/A[0,0]
    if phase_compensation ==  1:
        tTE3 = 1/A[0,0]*np.exp(-1j*k_0*(I.thickness_ML1)*np.tan(theta_out2)*np.sin(theta_out2))
        #/np.sqrt((Z_diec))
    else:
        tTE3 = 1/A[0,0]
    #(Z_diec2)
    A = multiLayerTransferMatrix(incidenceAngleRadians2, layerThickness, complexRelativePermittivity, frequency, 'tm')
    rTM3 = A[1,0]/A[0,0]
    tTM3 = 1/A[0,0]
    ############################################################################################################################

    return (tTE1* tTE2 * tTE3)




   
