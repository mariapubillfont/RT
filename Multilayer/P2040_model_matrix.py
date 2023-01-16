


import numpy as np
import matplotlib.pyplot as plt

    # Return final angle?
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
