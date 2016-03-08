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
