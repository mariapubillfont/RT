import numpy as np


def getReflectionCoefficients_multiLayer(k_0, layerThickness, polaritzation, complexPermittivity, incidentAngle):
    #k0 - free space wave number
    #layerthickness - thickness of each "slab"
    #polarization - TE or TM
    #complexPermittivity - vector of permittiviteies corresponding to each layer
    #incidentAngle - incident angle from air region

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
        A[i] = 0.5*np.exp(1j*k_n[i]*layerThickness[i]/np.cos(incidenceAngle_n[i]))*(A[i+1]*(1+Y)+B[i+1]*(1-Y))
        B[i] = 0.5*np.exp(-1j*k_n[i]*layerThickness[i]/np.cos(incidenceAngle_n[i]))*(A[i+1]*(1-Y)+B[i+1]*(1+Y))
        F[i] = 0.5*np.exp( 1j*k_n[i]*layerThickness[i]*np.cos(incidenceAngle_n[i]))*(F[i+1]*(1+W)+G[i+1]*(1-W))
        G[i] = 0.5*np.exp(-1j*k_n[i]*layerThickness[i]*np.cos(incidenceAngle_n[i]))*(F[i+1]*(1-W)+G[i+1]*(1+W))

    if polaritzation == 'TE':
        r = B[0]/A[0] #reflection and transmission for TE polaritzation
        t = 1/A[0]
    else:
        r = G[0]/F[0] #reflection and transmission for TM polartization
        t = 1/F[0]
    return t