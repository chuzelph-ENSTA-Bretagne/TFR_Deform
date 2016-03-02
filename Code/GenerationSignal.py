# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 19:10:39 2015

@author: Philippe
"""




import numpy as np
import math as m
import scipy
import pylab



def An(alpha,n):
    return 1j**n*np.exp(-1j*n*alpha)
    
# Cas du cylindre mou, champs nul aux niveau du contour en r
    
def Bn(alpha,k,r0,n):
    b=complex(scipy.special.jv(n, k*r0))
    a=complex(scipy.special.hankel1(n, k*r0))
    A=An(alpha,n)
    return -A*b/a
    
def GenSignal(k,r0,r,alpha,theta):
    # alpha : angle d'incidence de l'onde
    # theta : angle avec l'observateur
    # r     : distance de l'observateur au cylindre
    # r0    : rayon du cylindre
    # k     : w/c
    
    # Renvoie la liste des coefficient Bn de l'onde diffusée
    # Generation de Bn selon k*r0 et les composante du problème
    
    # Nombre de terme de la serie
    N = k*r0 + 10   

    # Calcule des coefficients     
    B=np.zeros(2*N+1,dtype=np.complex64)
    for j in range(len(B)):
        B[j]=Bn(alpha,k,r0,-N+j)
    return B
    
def Signal(k,r0,r,alpha,theta,tau,h):
    # Nombre de terme de la serie
    N = k*r0 + 50  
    t = scipy.linspace(0,tau,tau*h, endpoint=False)
    S = scipy.sin(250*2*m.pi*t)
    p = np.zeros(len(t),dtype=np.complex64)
    for j in range(len(t)):
        B = GenSignal(k,r0+r0*0.1*S[j],r,alpha,theta)
        tmp=complex(0)
        for l in range(len(B)):
            tmp=complex(tmp+B[l]*np.exp(1j*(-N+l)*theta)*scipy.special.hankel1(-N+l, k*r))
#        print '\n' + str(j) 
#        print tmp
        p[j]=tmp
    return p,t
        
def test ( a=5, b=10):
    print a+b
    
if __name__ == '__main__':
    fs = 800000.0        
    S,t = Signal(500,0.02,0.1,0,0,0.01,fs)
