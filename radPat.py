# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:29:39 2022

@author: mapf
"""
import numpy as np
import matplotlib.pyplot as plt
import input as I


#=============================================================================
def getUnitVector(x1,y1, x2, y2):
    vector=[x2-x1, y2-y1]
    norma = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    u = vector/norma
    return(u)
#=============================================================================


def norma(u):
    return np.sqrt(u[0]**2 + u[1]**2)

def getRadiationPattern(Ak_ap, path_length, nk, sk, dCk, Pk_ap):
    
    q = 0.1
    theta = np.linspace(0, np.pi, 400)
    R_obs =1.E9
    E = np.zeros(len(theta),complex)
    # print(Pk_ap)
    if 0:
        for j in range(0, len(theta)):
        
            xrj = np.cos(theta[j])*R_obs
            yrj = np.sin(theta[j])*R_obs
            
            E[j]=0.
            for k in range(0,len(Ak_ap)): 
                urkj_vector = getUnitVector(Pk_ap[k+1][0], Pk_ap[k+1][1], xrj, yrj)
                # plt.plot(Pk_ap[i+1][0], Pk_ap[i+1][1], 'x', color='green')
                rkj_distance = np.sqrt((xrj-Pk_ap[k+1,0])**2 + (yrj-Pk_ap[k+1,1])**2)
                # plt.quiver(Pk_ap[i+1][0], Pk_ap[i+1][1], rkj_vector[0], rkj_vector[1], color = 'blue')
                # plt.plot(Pk_ap[i+1][0], Pk_ap[i+1][1], 'x', color='blue')
                # print(rkj_distance)
                if 0:
                    fig = plt.figure(1)
                    fig.set_dpi(300)
                    plt.ylabel('rk')
                    plt.grid()
                
                Ek = sk[k+1,0]*urkj_vector[0] + sk[k+1,1]*urkj_vector[1]
                Ek = np.where(Ek>0,Ek**q,0.)
                
                E[j] += Ek*Ak_ap[k]*(np.exp(-1j*I.k0*(rkj_distance + path_length[k+1])))/rkj_distance\
                *(nk[k+1][0]*sk[k+1,0] + nk[k+1,1]*sk[k+1][1] - (nk[k+1,0]*urkj_vector[0]+nk[k+1,1]*urkj_vector[1]))*dCk[k]
                                                                                                                     
    if 1:
        Ez = np.zeros(len(theta)) + 0.j
        Xobs = np.cos(theta)*R_obs
        Yobs = np.sin(theta)*R_obs
        Xdip = Pk_ap[:,0]
        Ydip = Pk_ap[:,1]
        unvL = nk
        for ii in range(len(Pk_ap[:,0])):
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
            Ez[:] +=  Ak_ap[ii] * Fcos[:] * np.exp(-2j*np.pi*(path_length[ii]+dR[:])/I.wv) / (dR[:]) \
                    * (Fk + Gk[:]) * dCk[ii]
            # Ez[:] +=  Ak_ap[ii] * Fcos[:] * np.exp(-2j*np.pi*(path_length[ii]+dR[:])/I.wv) / (dR[:])    
        # return Ez                                                                                                                 ]
    
    return Ez, theta, dR