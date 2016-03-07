# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 19:10:39 2015

@author: Philippe
"""

import numpy as np
import math as m
import scipy as sc
import scipy.special as spe
import pylab as py


# Cas du cylindre mou, champs nul aux niveau du contour en r
    
def Bn(k,r0,n):
    b=complex(spe.jv(n, k*r0))
    a=complex(spe.hankel1(n, k*r0))
    A=1j**n
    return -A*b/a
    
def Field_scat(k,r0,r,alpha,theta,N):
    F=Bn(k,r0,0)*spe.hankel1(0, k*r)
    for n in range(1,N+1):
        F+=Bn(k,r0,n)*spe.hankel1(n, k*r)*2*np.cos(n*(theta-alpha))
    return F
    
def Signal(k,r,theta,r0,alpha,tau,npts,freq):
    # Variation du rayon du cylindre
    t = sc.linspace(0,tau,npts)
    tau_var=0.1 # strictement inférieur à 1
    V =r0*(1+tau_var*np.sin(2*m.pi*freq*t))
    
    # Nombre de terme de la serie
    N=int(np.ceil(k*r0*(1+tau_var)) + 10) 
    
    # Calcul champ incident
    Field_inc=np.exp(1j*k*r*np.cos(theta-alpha))
    
    # Calcul des champs diffusés
    Field=[Field_scat(k,r0var,r,alpha,theta,N) for r0var in V]+ Field_inc
    return Field,t
    


def Phase(rObjet,thetaObjet,alpha,r,theta):
    d1 = rObjet*np.cos(alpha-thetaObjet)
    rp=np.sqrt(rObjet**2+r**2-2*r*rObjet*np.cos(thetaObjet-theta))
    SthetaP = rObjet*np.sin(thetaObjet-theta)   # terme en 1/rp enleve
    CthetaP = (rp**2-r**2+rObjet**2)/(2*r)      # terme en 1/rp enleve
    thetaP=np.arctan2(SthetaP,CthetaP)
    return d1,rp,thetaP
    
def SignalObjDif(k,r,theta,r0,alpha,tau,npts,freq,rObjet,thetaObjet):
    d1,rp,thetaP=Phase(rObjet,thetaObjet,alpha,r,theta)
    S,t=Signal(k,rp,thetaP,r0,alpha,tau,npts,freq)
    return S*np.exp(1j*k*r*d1),t
    
    


    
if __name__ == '__main__':
    freq=250.
    tau = 5
    nbpts = 100*2*freq
    
    k=5
    r=20.
    theta=0.
    r0=1.
    alpha=0.
      
    # code de départ
    S,t = Signal(k,r,theta,r0,alpha,tau,nbpts,freq)
    
    py.plot(t,np.real(S))
    
    py.show()
    

    
    
    # code pour Chuzel Philippe, test le signal percu dans le cas ou le centre du cylindre n'est pas confondu avec l'origine du repere
    
    #d,r,the = Phase(2,2*np.pi/3,np.pi*(3/4),4,np.pi/4)
    #print(d,r,the)    
    #d,r,the = Phase(2,np.pi/4,np.pi*(3/4),4,np.pi/4)   
    #print('\n') 
    #print(d,r,the)    
    
    S,t = SignalObjDif(k,r,theta,r0,alpha,tau,nbpts,freq,4,np.pi/4)
    py.plot(t,np.real(S))
    
    py.show()
