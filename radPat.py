# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:29:39 2022

@author: mapf
"""
import numpy as np
import matplotlib.pyplot as plt
import input as I
import pandas as pd
from scipy.interpolate import interp1d

L = I.L
Array = I.Array
N = I.N
k0 = I.k0
alpha = I.alpha
beta = I.beta


# df = pd.read_excel('phaseDistribution_' + str(0) + '.xlsx', sheet_name='Sheet1')
# df_np = np.array(df)
# phase = df_np[:,1]
# x_phase = df_np[:,0]
# f = interp1d(x_phase, phase)
# xnew = (np.linspace(-L/2, L/2, num=N+1, endpoint=True))
# #phase_distribution = f(Array)    
# fig1 = plt.figure(1)
# plt.plot(x_phase, -phase, '.')
# fig1.set_dpi(300)
# plt.title('phase distribution/k0')
# plt.grid()
# plt.show() 


#=============================================================================
def getUnitVector(x1,y1, x2, y2):
    vector=[x2-x1, y2-y1]
    norma = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    u = vector/norma
    return(u)
#=============================================================================


def norma(u):
    return np.sqrt(u[0]**2 + u[1]**2)

def getRadiationPattern(Ak_ap, path_length, nk, sk, dCk, xap, yap, ts):
    q = 0.1
    theta = np.linspace(0, np.pi, 2500)
    R_obs =1.E9
    E = np.zeros(len(theta),complex)
    # print(Pk_ap)                                                                                                                    

    Ez = np.zeros(len(theta)) + 0.j
    Xobs = np.cos(theta)*R_obs
    Yobs = np.sin(theta)*R_obs
    Xdip = xap
    Ydip = yap
    unvL = nk
    for ii in range(len(xap)):
        dR = np.zeros(len(theta))
        vRx = Xobs[:]-Xdip[ii]
        vRy = Yobs[:]-Ydip[ii]
        dR = np.sqrt(vRx**2+vRy**2)  #Distance between sources and observation points.
        uvRx = vRx/dR; uvRy = vRy/dR # Unit vector associated with vR
        Gk = uvRx*unvL[ii,0] + uvRy*unvL[ii,1]
        uvSx, uvSy = sk[ii,0], sk[ii,1]
        Fk = uvSx*unvL[ii,0] + uvSy*unvL[ii,1]

        if 1: #feed == 'dipole':
            #uvSx, uvSy = np.cos(AoA[ii]), np.sin(AoA[ii])
            # uvSx, uvSy = uvs[0,ii], uvs[1,ii]
            # uvSx, uvSy = sk[ii,0], sk[ii,1]
            Fcos = uvSx*uvRx + uvSy*uvRy
            Fcos = np.where(Fcos<0,  0., Fcos**0.1)
        else:
            Fcos = np.ones_like(dR)
        if I.reflections == 1:
            Ez[:] +=  Ak_ap[ii] * Fcos[:] * np.exp(-1j*k0*(path_length[ii]+dR[:])) / (dR[:])*(Fk + Gk[:])*dCk[ii]*ts[ii]
        else:
            Ez[:] +=  Ak_ap[ii] * Fcos[:] * np.exp(-1j*k0*(path_length[ii]+dR[:])) / (dR[:])*(Fk + Gk[:])*dCk[ii]     
        # Ez[:] +=  Ak_ap[ii] * Fcos[:] * np.exp(-2j*np.pi*(path_length[ii]+dR[:])/I.wv) / (dR[:])    
    # return Ez                                                                                                                 ]
    
    return Ez, theta