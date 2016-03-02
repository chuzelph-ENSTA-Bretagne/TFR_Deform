# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 09:36:15 2015

@author: Philippe
"""

from scipy.fftpack import fft, fftfreq
from scipy.signal import chirp
import matplotlib.pyplot as plt   
import numpy as np
import math as m
import scipy

'''
Bibliotheques personnalisées
'''   

import GenerationSignalV2
import myspectrogram2
#import WiegnerVille
import CodeFournis

'''
Début des fonction
'''

'''
génère les coefficient An, vrai dans tous les cas
'''

def An(alpha,n):
    return 1j**n*np.exp(-1j*n*alpha)

'''
Définition des Bn    OK A NE PAS TOUCHER!!!!!
'''
    
# Cas du cylindre mou, champs nul aux niveau du contour en r0
def Bn(alpha,k,r0,n):
    b=complex(scipy.special.jv(n, k*r0))
    a=complex(scipy.special.hankel1(n, k*r0))
    A=An(alpha,n)
    return -A*b/a
    
    
''' 
Definition des des Bn a un t donné    OK A NE PAS TOUCHER!!!!!
'''

def GenSignal(k,r0,alpha,theta,N):
    # alpha : angle d'incidence de l'onde
    # theta : angle avec l'observateur
    # r     : distance de l'observateur au cylindre
    # r0    : rayon du cylindre
    # k     : w/c
    # N     : nombre de terme de la serie
    
    # Renvoie la liste des coefficient Bn de l'onde diffusée
    # Generation de Bn selon k*r0 et les composante du problème
    
    # Calcule des coefficients     
    B=np.zeros(2*N+1,dtype=np.complex64)
    for j in range(len(B)):
        B[j]=Bn(alpha,k,r0,-N+j)
    return B
    
''' 
Definition du signal    OK A NE PAS TOUCHER!!!!!
'''

def Signal(k,r0,r,alpha,theta,tau,h):
    # Nombre de terme de la serie
    N = k*r0 + 30  
    t = scipy.linspace(0,tau,tau*h, endpoint=False)
    S = scipy.sin(250*2*m.pi*t)
    p = np.zeros(len(t),dtype=np.complex64)
    for j in range(len(t)):
        B = GenSignal(k,r0+r0*0.1*S[j],alpha,theta,N)
        tmp=complex(0)
        for l in range(len(B)):
            tmp=complex(tmp+B[l]*np.exp(1j*(-N+l)*theta)*scipy.special.hankel1(-N+l, k*r))
        p[j]=tmp
    return p,t
    
def SignalEllipse(k,a,b,w0,r,alpha,theta,tau,h):
    # Nombre de terme de la serie
    r0=m.sqrt((a*a+b*b)/2)
    N = k*r0 + 30  
    
    t = scipy.linspace(0,tau,tau*h, endpoint=False)
    p = np.zeros(len(t),dtype=np.complex64)
    for j in range(len(t)):
        #Calcul du r(t)
        rt=m.sqrt(a*a*m.cos(2*m.pi*t[j]*w0+alpha)**2+b*b*m.sin(2*m.pi*t[j]*w0+alpha)**2)
        B = GenSignal(k,rt,alpha,theta,N)
        tmp=complex(0)
        for l in range(len(B)):
            tmp=complex(tmp+B[l]*np.exp(1j*(-N+l)*theta)*scipy.special.hankel1(-N+l, k*r))
        p[j]=tmp
    return p,t
    
def GenerationSignalSimple(tau,h,w1,w2,t1,t2):
    t = np.arange(0.0, tau, h)
    s1 = 10*scipy.sin(2*np.pi*w1*t)
    s2 = 10*scipy.sin(2*np.pi*w2*t)
    
    # create a transient "chirp"
    mask = np.where(np.logical_and(t>t1, t<t2), 1.0, 0.0)
    s2 = s2 * mask
    
    # add some noise into the mix
    nse = 0.01*np.random.randn(len(t))
    x = s1 + s2 + nse # the signal
    return x

'''
Fonctions d'affichage    OK A NE PAS TOUCHER!!!!!
'''
    
def AffichageFFT(x,fs,tau):
    FenAcq = x.size             # taille de la fenetre temporelle
    t=scipy.linspace(0, tau, x.size)
    
    # calcul de la TFD par l'algo de FFT
    signal_FFT = abs(fft(x))    # on ne recupere que les composantes reelles
    
    # recuperation du domaine fréquentiel
    signal_freq = fftfreq(FenAcq,1/float(fs))
    
    # extraction des valeurs reelles de la FFT et du domaine frequentiel
    signal_FFT = signal_FFT[0:len(signal_FFT)//2]
    signal_freq = signal_freq[0:len(signal_freq)//2]  
    
    #affichage du signal
    plt.figure()
    plt.subplot(211)
    plt.title('Signal et son spectre')
    plt.plot(t, x)
    plt.xlabel('Temps (s)')
    plt.ylabel('Amplitude')
    plt.show()
    
    #affichage du spectre du signal
    plt.subplot(212)
    plt.plot(signal_freq,signal_FFT)
    plt.xlabel('Frequence (Hz)') 
    plt.ylabel('Amplitude')
    
    
def genSpectrogram(x,fs,wlen=1024,h=8,nfft=1024):
    X , T , F = myspectrogram2.myspectrogram(x, wlen, h, nfft, fs)
    plt.figure()
    plt.imshow(abs(X), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
    plt.xlabel('Time')
    plt.ylabel('Frequency')
    plt.colorbar()
    yl= np.float32(np.linspace(0, len(F)-1, 11))
    ylocs = np.float32(np.round(np.linspace(0, F[len(F)-1], 11)))
    plt.yticks(yl, ["%.02f" % i for i in ylocs])
    

    plt.show()
    return T,F,X
 
 
'''
fonctions données --- Partie à supprimer!
'''


    
'''
#######################  MAIN  ############################
'''


if __name__ == '__main__':
       
    fe=500
    tau=1
    t1=0.1
    t2=0.5
    f1=200
    f2=150
    x=GenerationSignalSimple(tau,1/float(fe),f1,f2,t1,t2)
    AffichageFFT(x,fe,tau)
    X , T , F = genSpectrogram(x,fe,wlen=32,h=8,nfft=1024)
    

#    fe=500
#    fr=50
#    tau = 5
#    nbpts = tau*fe
#    k=5
#    r=20
#    theta=m.pi #m.pi/4
#    r0=1
#    alpha=0
#    x,t = GenerationSignalV2.Signal(k,r,theta,r0,alpha,tau,nbpts,fr) 
#    b=np.random.randn(len(x))*0.01
#    x=abs(x)-np.mean(abs(x))
#    AffichageFFT(x,fe,tau)
#    T2 , F2 , X2 = genSpectrogram(x,fe,wlen=32,h=1,nfft=1024)
    
    '''
    Partie test de fonction avec signal chirp
    '''
    
#    fe = 1000
#    tau = 0.5
#    t = scipy.linspace(0,tau,fe*tau)
#    x = chirp(t, 0, tau, 300) + chirp(t, 350, tau, 50)
#    AffichageFFT(x,fe,tau)
#    T2 , F2 , X2 = genSpectrogram(x,fe,wlen=32,h=1,nfft=1024)
    
    '''
    test de Wiegner-Ville
    '''
    
#    X = WiegnerVille.DWVD(x,fe,tau)
#    plt.figure()
#    plt.imshow(abs(X), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
#    plt.xlabel('Time')
#    plt.ylabel('Frequency')
#    plt.colorbar()
#    plt.show()

    '''
    partie test du code fournie
    '''
    
#    test = CodeFournis.stft(x, 1024)
#    plt.figure()
#    plt.imshow(np.transpose(abs(test)), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
#    plt.xlabel('Time')
#    plt.ylabel('Frequency')
#    plt.colorbar()
    
    # X2 , T2 , F2 = CodeFournis.wvd(np.transpose(x),N=1024)
    # plt.figure()
    # plt.imshow(np.transpose(abs(X2)), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
    # plt.xlabel('Time')
    # plt.ylabel('Frequency')
    # plt.colorbar()
    # plt.show()
    
    '''
    partie poubelle
    '''
    # X3 , T3 , F3 = CodeFournis.tfrspwv(np.transpose(x),N=1024)
    # plt.figure()
    # plt.imshow(np.transpose(abs(X3)), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
    # plt.xlabel('Time')
    # plt.ylabel('Frequency')
    # plt.colorbar()
    # plt.show()
    
