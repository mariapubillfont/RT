
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
