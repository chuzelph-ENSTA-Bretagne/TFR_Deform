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
from scipy.fftpack import fft, fftfreq
from scipy.signal import chirp
import myspectrogram
import tftb

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
    tau_var=0.2 # strictement inférieur à 1
    V =r0*(1+tau_var*np.sin(2*m.pi*freq*t))
    
    # Nombre de terme de la serie
    N=int(np.ceil(k*r0*(1+tau_var)) + 10) 
    
    # Calcul champ incident
    Field_inc=np.exp(1j*k*r*np.cos(theta-alpha))
    
    # Calcul des champs diffusés
    Field=[Field_scat(k,r0var,r,alpha,theta,N) for r0var in V]+ Field_inc
    py.figure()
    py.plot(t,V)
    py.title('rayon en cours du temps')
    py.figure()
    py.plot(t,np.real(Field-Field_inc))
    py.title('champs diffuse')
    py.figure()
    py.plot(t,np.real(np.ones((npts,1))*Field_inc))
    py.title('champs incident')
    return Field,t


def AffichageFFT(x,fs,tau):
    FenAcq = x.size             # taille de la fenetre temporelle
    t=sc.linspace(0, tau, x.size)
    
    # calcul de la TFD par l'algo de FFT
    signal_FFT = np.real(fft(x))    # on ne recupere que les composantes reelles de x
    
    # recuperation du domaine fréquentiel
    signal_freq = fftfreq(FenAcq,1/float(fs))
    
    # extraction des valeurs reelles de la FFT et du domaine frequentiel
    signal_FFT = signal_FFT[0:len(signal_FFT)//2]
    signal_freq = signal_freq[0:len(signal_freq)//2]  
    
    #affichage du signal
    py.figure()
    py.subplot(211)
    py.title('Signal et son spectre')
    py.plot(t, x)
    py.xlabel('Temps (s)')
    py.ylabel('Amplitude')
    
    #affichage du spectre du signal
    py.subplot(212)
    py.plot(signal_freq,signal_FFT)
    py.xlabel('Frequence (Hz)') 
    py.ylabel('Amplitude')
    py.show()
    
    
def genSpectrogram(x,fs,wlen=1024,h=8,nfft=1024):
    X , T , F = myspectrogram.myspectrogram(x, wlen, h, nfft, fs)
    py.figure()
    py.imshow(abs(X), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
    py.xlabel('Time')
    py.ylabel('Frequency')
    py.colorbar()
    py.show()
    return T,F,X

    
if __name__ == '__main__':
    freq=500.
    tau = 0.05
    nbpts = 20000*tau*freq
    
    k=100
    r=20.
    theta=0.
    r0=1.
    alpha=0.
      
    # code de départ
    S,t = Signal(k,r,theta,r0,alpha,tau,nbpts,freq)
    S = np.array(S)
    AffichageFFT(np.real(S-np.mean(np.real(S))),20*freq,tau)
    #T2 , F2 , X2 = genSpectrogram(np.real(S-np.mean(S)),200*freq,wlen=2048,h=50,nfft=2048)
    tfr = tftb.processing.ShortTimeFourierTransform(np.real(S-np.mean(S)),n_fbins = 2048,fwindow = np.hamming(512))
    tfr.run()
    tfr.plot(kind='cmap', show_tf=True)
