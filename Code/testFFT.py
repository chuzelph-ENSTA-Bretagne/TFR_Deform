# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 18:24:15 2015

@author: Philippe
"""

from numpy import pi, sin, linspace, log10, random
from scipy.fftpack import fft, fftfreq
import matplotlib.pyplot as plt      
import myspectrogram2
#import pylab
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
# définition des constantes du signal
K = 2*pi       # facteur conversion période/fréquence
A0 = 4         # amplitude fréquence fondamentale
A1 = 8         # amplitude première harmonique
f0 = 2         # fréquence fondamentale (Hz)
f1 = 8         # fréquence première harmonique (Hz)

# définition temporelle du signal
t0 = 0         # début de l'acquisition du signal
t1 = 10        # fin de l'acquisition (s)

# définition des paramètres d'échantillonnage
FreqEch = 1024                  # fréquence d'échantillonage
PerEch = 1./FreqEch             # période d'échantillonnage
N = FreqEch*(t1 - t0)           # nombre de points échantillonnés sur l'intervalle

# définition du temps
t = linspace(t0, t1, N)

# définition du signal
signal = A0*sin(f0*K*t) + A1*sin(f1*K*t)

# bruitage aléatoire du signal composite. Le signal composite est
# noyé dans le bruit.
Ab = 20                         # amplitude du bruit
fb = 50                         # fréquence du bruit
#signal = signal + Ab*random.random(len(signal))*sin(fb*K*t) 
signal = signal + Ab*sin(random.random(len(signal))*fb*K*t)
# définition des données de FFT
FenAcq = signal.size             # taille de la fenetre temporelle
    
# calcul de la TFD par l'algo de FFT
signal_FFT = abs(fft(signal))    # on ne récupère que les composantes réelles

# récupération du domaine fréquentiel
signal_freq = fftfreq(FenAcq,PerEch)

# extraction des valeurs réelles de la FFT et du domaine fréquentiel
signal_FFT = signal_FFT[0:len(signal_FFT)//2]
signal_freq = signal_freq[0:len(signal_freq)//2]

#affichage du signal
plt.figure(1)
plt.subplot(211)
plt.title('Signal et son spectre')
plt.ylim(-(Ab+10), Ab+10)
plt.plot(t, signal)
plt.xlabel('Temps (s)'); plt.ylabel('Amplitude')
plt.show()

#affichage du spectre du signal
plt.subplot(212)
plt.xlim(0,fb+5)
plt.plot(signal_freq,signal_FFT)
plt.xlabel('Frequence (Hz)'); plt.ylabel('Amplitude')


dt = 0.0005
t = np.arange(0.0, 20.0, dt)
s1 = 10*sin(2*pi*100*t)
s2 = 10*sin(2*pi*400*t)

# create a transient "chirp"
mask = np.where(np.logical_and(t>10, t<12), 1.0, 0.0)
s2 = s2 * mask

# add some noise into the mix
nse = 0.01*np.random.randn(len(t))

x = s1 + s2 + nse # the signal

X , T , F = myspectrogram2.myspectrogram(x, wlen=1024, h=8, nfft=1024, fs=1/dt)
print X.shape
print len(T)
print len(F)
plt.figure(2)
plt.imshow(abs(X), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
plt.xlabel('Time')
plt.ylabel('Frequency')
plt.show()

#fig = plt.figure(3)
#ax = plt.axes(projection='3d')
#ax.plot_surface(np.outer(F, np.ones(len(T))),np.outer(T, np.ones(len(F))).T,abs(X), cmap='jet')
#plt.show()

