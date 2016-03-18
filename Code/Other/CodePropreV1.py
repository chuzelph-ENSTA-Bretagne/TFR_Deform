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

import GenerationSignal
import myspectrogram

'''
Début des fonction
'''

    
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
    
    #affichage du spectre du signal
    plt.subplot(212)
    plt.plot(signal_freq,signal_FFT)
    plt.xlabel('Frequence (Hz)') 
    plt.ylabel('Amplitude')
    plt.show()
    
    
def genSpectrogram(x,fs,wlen=1024,h=8,nfft=1024):
    X , T , F = myspectrogram.myspectrogram(x, wlen, h, nfft, fs)
    plt.figure()
    plt.imshow(abs(X), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
    plt.xlabel('Time')
    plt.ylabel('Frequency')
    plt.colorbar()
    plt.show()
    return T,F,X
 

    
'''
#######################  Test unitaire  ############################
'''


if __name__ == '__main__':
       
    fe=500
    tau=1
    t1=0.1
    t2=0.5
    f1=200
    f2=150
    #x=GenerationSignalSimple(tau,1/float(fe),f1,f2,t1,t2)
    #AffichageFFT(x,fe,tau)
    #X , T , F = genSpectrogram(x,fe,wlen=32,h=8,nfft=1024)
    

    fe=1000
    fr=50
    tau = 1
    nbpts = tau*fe
    k=5
    r=20
    theta=m.pi #m.pi/4
    r0=1
    alpha=0
    x,t = GenerationSignal.Signal(k,r,theta,r0,alpha,tau,nbpts,fr) 
    x=abs(x)-np.mean(abs(x))
    AffichageFFT(x,fe,tau)
    T2 , F2 , X2 = genSpectrogram(x,fe,wlen=32,h=1,nfft=1024)
    
    '''
    Partie test de fonction avec signal chirp
    '''
    
#    fe = 1000
#    tau = 0.5
#    t = scipy.linspace(0,tau,fe*tau)
#    x = chirp(t, 0, tau, 300) + chirp(t, 350, tau, 50)
#    AffichageFFT(x,fe,tau)
#    T2 , F2 , X2 = genSpectrogram(x,fe,wlen=32,h=1,nfft=1024)

   
    
