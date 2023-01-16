
import numpy as np
import matplotlib.pyplot as plt
# f = 13
# c0 = 0.299792458 
# wv = c0/f
# k_0 = 2*np.pi/(c0/f)
#a_in = [0 ,40]
a_in = np.deg2rad(np.arange(0, 89, 1))

# General parameters
frequency = 3e8
lambd = 3e8/frequency
wv = lambd
k_0 = 2*np.pi/lambd
d_core = lambd # Core thickness
#thickness = np.arange(0,0.1,0.00001)
#n2 = 5.31
c = 0.0326
d = 0.8095
cond = c*frequency**d #S/m
er = 2.5 -1j*10
#permittivity = n2 - 1j*(17.98*cond/f)
#tan_delta = 0.00066
tan_delta = 10
permittivity = er*(1-1j*tan_delta)
thickness = d_core


TML = np.zeros(len(a_in), dtype=np.complex_)
RML = np.zeros(len(a_in),dtype=np.complex_)

def getReflectionCoefficients(wv, thickness, polaritzation, permittivity, incidentAngle):
    eta = np.sqrt(permittivity-np.sin(incidentAngle)**2)
    if polaritzation == 'TE': R_aux = (np.cos(incidentAngle)-eta)/(np.cos(incidentAngle)+eta)
    else: R_aux = (permittivity*np.cos(a_in)-eta)/(permittivity*np.cos(incidentAngle)+eta)
    q = eta*2*np.pi*thickness/wv
    T = ((1-R_aux**2)*np.exp(-1j*q))/(1-R_aux**2*np.exp(-2j*q))
    R = R_aux*(1-np.exp(-2j*q))/(1-R_aux**2*np.exp(-2j*q))
    return T, R

layerThickness = [0, wv/(np.sqrt(np.sqrt(er))*4), d_core, wv/(np.sqrt(np.sqrt(er))*4), 0]
complexPermittivity = [1, np.sqrt(er), er, np.sqrt(er), 1]

# layerThickness = [0, wv/(np.sqrt(np.sqrt(er))*4), 0]
# complexPermittivity = [1, np.sqrt(er), er]

# layerThickness = [0, thickness, 0]
# complexPermittivity = [1, er, 1]

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
        A =  np.ones(N, dtype=np.complex_)/np.sqrt(complexPermittivity[N-1])
        B = np.zeros(N, dtype=np.complex_)
        F =  np.ones(N, dtype=np.complex_)/np.sqrt(complexPermittivity[N-1])
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


TML, RML = getReflectionCoefficients_multiLayer(k_0, layerThickness, 'TE', complexPermittivity, a_in)


# T, R = getReflectionCoefficients(wv, thickness, 'TM', permittivity, a_in)

# plt.plot(thickness, 20*np.log10(abs(R)))
# plt.plot(thickness, 20*np.log10(abs(T)))
# plt.xlabel('Slab thickness (m)')
# print(20*np.log10(abs(R)), 20*np.log10(abs(T)))
# # plt.plot(c0/wv_range, 20*np.log10(abs(R)))
# # plt.plot(c0/wv_range, 20*np.log10(abs(T)))
# #plt.xlabel('Frequency (GHz)')
# plt.ylabel('T and R coefficents amplitude (dB)')

# plt.title('Transmission and reflection coefficents')
# plt.legend(['Reflection', 'Tramission'])

# plt.grid()
# plt.show()


plt.figure()
plt.ylim([-20, 1.1])
plt.plot(np.rad2deg(a_in),  20*np.log10(abs(RML)))
plt.plot(np.rad2deg(a_in),  20*np.log10(abs(TML)))
plt.plot(np.rad2deg(a_in), 10*np.log10(abs(TML)**2 + abs(RML)**2))
plt.legend(['Reflection', 'Transmission'])
plt.grid()
plt.show()



# plt.figure()
# plt.plot(np.rad2deg(a_in), abs(R))
# plt.plot(np.rad2deg(a_in), abs(T))
# plt.plot(np.rad2deg(a_in), abs(R)**2 + abs(T)**2)
# plt.legend(['Reflection', 'Tramission'])
# plt.show()

# # % [rTE,rTM,tTE,tTM] = P2040multilayer(k_0,a_in,layerThickness,complexRelativePermittivity)
# # % This function calculates the transmission and reflection of a plane wave through single- or multi-layer materials
# # % according to the model in Rec. ITU-R P.2040 Effects of building materials and structures on radiowave propagation
# # % above about 100 MHz available at http://www.itu.int/rec/R-REC-P.2040-1-201507-I.
# # %
# # % The input parameters are:
# # %    k_0                              A column vector with the free space wave numbers
# # %    incidenceAngle                   The incidence angle in radians with respect to the normal
# # %    layerThickness                   A row vector with the layer a_incknesses in m
# # %    complexRelativePermittivity      A matrix with the complex relative permittivity per layer and frequency
# # %
# # % Output parameters are the reflection coefficients and the transmission coefficients for the TE and TM modes
# f = 13e9
# c = 3e8
# k_0 = 2*np.pi/(c/f)
# a_in = [0, 0.6981]
# layerThickness = [0.1]
# complexRelativePermittivity = 2.5
# #def P2040multilayer(k_0,a_in,layerThickness,complexRelativePermittivity): 
# if len(layerThickness)==1:
#     #Single material, use eqs (37ab), (43ab), and (44) from P.2040
#     cosTheta = np.cos(a_in)
#     sinTheta = np.sin(a_in)
#     eta = complexRelativePermittivity
    
#     sqrtEtaSin = np.sqrt(eta-sinTheta**2)
#     ReTE =      (cosTheta-sqrtEtaSin)/(cosTheta+sqrtEtaSin)
#     ReTM = (eta*cosTheta-sqrtEtaSin)/(eta*cosTheta+sqrtEtaSin)
#     # % Note: these two are not needed here
#     # % TeTE =                 2*cosTheta./(cosTheta+sqrtEtaSin);
#     # % TeTM =      2*sqrt(eta).*cosTheta./(eta.*cosTheta+sqrtEtaSin);
    
#     q = layerThickness*sqrtEtaSin*k_0
#     rTE =    ReTE*(1-np.exp(-2j*q))/(1-ReTE**2*np.exp(-2j*q))
#     rTM =    ReTM*(1-np.exp(-2j*q))/(1-ReTM**2*np.exp(-2j*q))
#     tTE = (1-ReTE**2)*np.exp(-1j*q)/(1-ReTE**2*np.exp(-2j*q))
#     tTM = (1-ReTM**2)*np.exp(-1j*q)/(1-ReTM**2*np.exp(-2j*q))
    
# else:
#     # % Multi-layer material, use eqs (39)-(42) from P.2040
    
#     # #% FIXME: For efficiency, the following two checks could be removed if the air layers are always (or never) present
#     # if any(complexRelativePermittivity(:,1)~=1)
#     #     % Add a zero thickness air layer in front
#     #     layerThickness = [0 layerThickness];
#     #     complexRelativePermittivity = [ones(length(k_0),1) complexRelativePermittivity];
#     # end
#     # if any(complexRelativePermittivity(:,end)~=1)
#     #     % Add a zero thickness air layer in at back
#     #     layerThickness = [layerThickness 0];
#     #     complexRelativePermittivity = [complexRelativePermittivity ones(length(k_0),1)];
#     # end
    
#     nLayers = len(layerThickness)
    
#     eta_n = complexRelativePermittivity
#     k_n = k_0*np.sqrt(eta_n)
    
#     incidenceAngle_n = np.asin(np.sin(a_in)/np.sqrt(eta_n))
    
#     N = nLayers
    
#     # % Initialize
#     A =  np.ones(len(k_0),1)
#     B = np.zeros(len(k_0),1)
#     F =  np.ones(len(k_0),1)
#     G = np.zeros(len(k_0),1)
    
#     # % Calculate backwards from last layer to first
#     #for N = (nLayers-1):-1:1
#     # for int i in range(N, 0):
#     #     W = cos(incidenceAngle_n(:,N+1))./cos(incidenceAngle_n(:,N)).*sqrt(eta_n(:,N)./eta_n(:,N+1));
#     #     Y = cos(incidenceAngle_n(:,N+1))./cos(incidenceAngle_n(:,N)).*sqrt(eta_n(:,N+1)./eta_n(:,N));
#     #     A(:,N) = 0.5*exp( 1i*k_n(:,N).*layerThickness(:,N).*cos(incidenceAngle_n(:,N))).*(A(:,N+1).*(1+Y)+B(:,N+1).*(1-Y));
#     #     B(:,N) = 0.5*exp(-1i*k_n(:,N).*layerThickness(:,N).*cos(incidenceAngle_n(:,N))).*(A(:,N+1).*(1-Y)+B(:,N+1).*(1+Y));
#     #     F(:,N) = 0.5*exp( 1i*k_n(:,N).*layerThickness(:,N).*cos(incidenceAngle_n(:,N))).*(F(:,N+1).*(1+W)+G(:,N+1).*(1-W));
#     #     G(:,N) = 0.5*exp(-1i*k_n(:,N).*layerThickness(:,N).*cos(incidenceAngle_n(:,N))).*(F(:,N+1).*(1-W)+G(:,N+1).*(1+W));
#     # end
#     # rTE = B(:,1)./A(:,1);
#     # rTM = G(:,1)./F(:,1);
#     # tTE = 1./A(:,1);
#     # tTM = 1./F(:,1);
# #return 0

# #return rTE,rTM,tTE,tTM
