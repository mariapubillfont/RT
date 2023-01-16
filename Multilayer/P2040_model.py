
import numpy as np
import matplotlib.pyplot as plt


def P2040multilayer(k_0,incidenceAngleRadians,layerThickness,complexRelativePermittivity):
# [rTE,rTM,tTE,tTM] = P2040multilayer(k_0,incidenceAngleRadians,layerThickness,complexRelativePermittivity)
# This function calculates the transmission and reflection of a plane wave through single- or multi-layer materials
# according to the model in Rec. ITU-R P.2040 Effects of building materials and structures on radiowave propagation
# above about 100 MHz available at https://www.itu.int/rec/R-REC-P.2040-2-202109-I.
#
# The input parameters are:
#    k_0                              A column vector with the free space wave numbers
#    incidenceAngle                   The incidence angle in radians with respect to the normal
#    layerThickness                   A row vector with the layer thicknesses in m
#    complexRelativePermittivity      A matrix with the complex relative permittivity per layer and frequency
#
# Output parameters are the reflection coefficients and the transmission coefficients for the TE and TM modes



    if np.size(layerThickness)==1:
        # Single material, use eqs (37ab), (43ab), and (44) from P.2040
        cosTheta = np.cos(incidenceAngleRadians)
        sinTheta = np.sin(incidenceAngleRadians)
        eta = complexRelativePermittivity
        
        np.sqrtEtaSin = np.sqrt(eta-sinTheta^2)
        ReTE =      (cosTheta-np.sqrtEtaSin)/(cosTheta+np.sqrtEtaSin)
        ReTM = (eta*cosTheta-np.sqrtEtaSin)/(eta*cosTheta+np.sqrtEtaSin)
        # Note: these two are not needed here
        # TeTE =                 2*cosTheta/(cosTheta+np.sqrtEtaSin)
        # TeTM =      2*np.sqrt(eta)*cosTheta/(eta*cosTheta+np.sqrtEtaSin)
        
        q = k_0*layerThickness*np.sqrtEtaSin
        rTE =    ReTE*(1-np.exp(-2j*q))/(1-ReTE^2*np.exp(-2j*q))
        rTM =    ReTM*(1-np.exp(-2j*q))/(1-ReTM^2*np.exp(-2j*q))
        tTE = (1-ReTE^2)*np.exp(-1j*q)/(1-ReTE^2*np.exp(-2j*q))
        tTM = (1-ReTM^2)*np.exp(-1j*q)/(1-ReTM^2*np.exp(-2j*q))
        
    else:
        # Multi-layer material, use eqs (39)-(42) from P.2040
        
    #     # FIXME: For efficiency, the following two checks could be removed if the air layers are always (or never) present
    #     if any(complexRelativePermittivity(:,1)~=1)
    #         # Add a zero thickness air layer in front
    #         layerThickness = [0 layerThickness]
    #         complexRelativePermittivity = [ones(length(k_0),1) complexRelativePermittivity]
    #     end
    #     if any(complexRelativePermittivity(:,end)~=1)
    #         # Add a zero thickness air layer in at back
    #         layerThickness = [layerThickness 0]
    #         complexRelativePermittivity = [complexRelativePermittivity ones(length(k_0),1)]
    #     end
        
        nLayers = np.size(layerThickness)
        
        eta_n = complexRelativePermittivity
        k_n = k_0*np.sqrt(eta_n)
        
        incidenceAngle_n = np.arcsin(np.sin(incidenceAngleRadians)/np.sqrt(eta_n))
        
        # Initialize
        A =  np.ones([np.size(k_0),nLayers], dtype=np.complex_) / (np.cos(incidenceAngle_n[:,nLayers-1]) * np.sqrt(eta_n[:,nLayers-1]))
        B = np.zeros([np.size(k_0),nLayers], dtype=np.complex_)
        F =  np.ones([np.size(k_0),nLayers], dtype=np.complex_) / (np.cos(incidenceAngle_n[:,nLayers-1]) / np.sqrt(eta_n[:,nLayers-1]))
        G = np.zeros([np.size(k_0),nLayers], dtype=np.complex_)
        
        # Calculate backwards from last layer to first
        for N in range (nLayers-1,-1,-1):
            Y = np.cos(incidenceAngle_n[:,N+1])/np.cos(incidenceAngle_n[:,N])*np.sqrt(eta_n[:,N+1]/eta_n[:,N])
            A[:,N] = 0.5*np.exp( 1j*k_n[:,N]*layerThickness[:,N]*np.cos(incidenceAngle_n[:,N]))*(A[:,N+1]*(1+Y)+B[:,N+1]*(1-Y))
            B[:,N] = 0.5*np.exp(-1j*k_n[:,N]*layerThickness[:,N]*np.cos(incidenceAngle_n[:,N]))*(A[:,N+1]*(1-Y)+B[:,N+1]*(1+Y))
            W = np.cos(incidenceAngle_n[:,N+1])/np.cos(incidenceAngle_n[:,N])*np.sqrt(eta_n[:,N]/eta_n[:,N+1])
            F[:,N] = 0.5*np.exp( 1j*k_n[:,N]*layerThickness[:,N]*np.cos(incidenceAngle_n[:,N]))*(F[:,N+1]*(1+W)+G[:,N+1]*(1-W))
            G[:,N] = 0.5*np.exp(-1j*k_n[:,N]*layerThickness[:,N]*np.cos(incidenceAngle_n[:,N]))*(F[:,N+1]*(1-W)+G[:,N+1]*(1+W))
        
        rTE = B[:,1]/A[:,1]
        rTM = G[:,1]/F[:,1]
        tTE = 1/(np.cos(incidenceAngle_n[:,nLayers])/np.cos(incidenceAngle_n[:,1]))/A[:,1]*np.sqrt(eta_n[:,nLayers]/eta_n[:,1])
        tTM = 1/(np.cos(incidenceAngle_n[:,nLayers])/np.cos(incidenceAngle_n[:,1]))/F[:,1]*np.sqrt(eta_n[:,nLayers]/eta_n[:,1])
    return rTE,rTM,tTE,tTM