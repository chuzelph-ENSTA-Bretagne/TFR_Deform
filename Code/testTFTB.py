# -*- coding: utf-8 -*-


from scipy.fftpack import fft, fftfreq
from scipy.signal import chirp, hilbert,find_peaks_cwt
import matplotlib.pyplot as plt   
import numpy as np
import math as m
import scipy
import tftb
try:
    from skimage import filters
except ImportError:
    from skimage import filter as filters
from skimage import exposure

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

'''
Fonctions d'affichage    OK A NE PAS TOUCHER!!!!!
'''
    
def AffichageFFT(x,fs,tau):
    FenAcq = x.size             # Nombre de point pour la FFT
    t=scipy.linspace(0, tau, x.size)
    
    # calcul de la TFD par l'algo de FFT
    signal_FFT = abs(fft(x))    # on ne recupere que les composantes reelles
    
    # recuperation du domaine fréquentiel
    signal_freq = fftfreq(FenAcq,1/float(fs))
    
    # extraction des valeurs reelles de la FFT et du domaine frequentiel
    signal_FFT = signal_FFT[0:len(signal_FFT)//2]
    signal_freq = signal_freq[0:len(signal_freq)//2]  
    
    #affichage du signal
    plt.figure('Affiche du signal temporelle et de sa densite spectrale')
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

''' 
Definition d'un signal de test   OK A NE PAS TOUCHER!!!!!
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
    return x,t

def lawFrequency(tfr,t,f):
	tmp = np.zeros((len(f),len(t)))
	k=0
	print(np.shape(tfr))
	print(np.shape(tmp))
	print(np.shape(t))
	print(np.shape(f))

	# création du chemin
	for i in range(len(t)):
		if(i==0):
			peakind = find_peaks_cwt(tfr[:,0], np.arange(1,10))
			print(peakind[0])
			tmp[peakind[0],0] = 1
			k = peakind[0]
			print(k)
		else:
			Vartmp = 0
			Indice = 0
			for j in [-1,0,1]:
				if(Vartmp<tmp[k+j,i-1]):
					Vartmp = tmp[k+j,i]
					Indice = j
			k = k + Indice
			print(k)
			tmp[k,i] = 1

	plt.figure()
	plt.imshow(tmp, origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
	plt.xlabel('Time')
	plt.ylabel('Frequency')
	plt.colorbar()
	plt.show()


	#print(peakind, freqs[peakind][0], tfr.tfr[peakind,30][0])

if __name__ == '__main__':

	fe=500
	tau=1
	t1=0.2
	t2=0.7
	f1=200
	f2=150
	x , t =GenerationSignalSimple(tau,1/float(fe),f1,f2,t1,t2)
	AffichageFFT(x,fe,tau)

	tfr = tftb.processing.ShortTimeFourierTransform(x)
	tfr.run()
	tfr.plot(kind='cmap', show_tf=True)
	freqs = tfr.freqs
	

	lawFrequency(tfr.tfr,tfr.ts,freqs)

	val = filters.threshold_otsu(tfr.tfr)

	hist, bins_center = exposure.histogram(tfr.tfr)
	print(val)
	plt.figure(figsize=(9, 4))
	plt.subplot(131)
	plt.imshow(tfr.tfr, cmap='gray', interpolation='nearest')
	plt.axis('off')
	plt.subplot(132)
	plt.imshow(tfr.tfr < val, cmap='gray', interpolation='nearest')
	plt.axis('off')
	plt.subplot(133)
	plt.plot(bins_center, hist, lw=2)
	plt.axvline(val, color='k', ls='--')

	plt.tight_layout()
	plt.show()

	#wvd = tftb.processing.WignerVilleDistribution(x)
	#wvd.run()
	#wvd.plot(kind='cmap', show_tf=True)   #kind='contour', show_tf=True
    
	#wvd2 = tftb.processing.cohen.PseudoWignerVilleDistribution(hilbert(x))
	#wvd2.run()
	#wvd2.plot(kind='contour',show_tf=True)   #kind='contour', show_tf=True


	#wvd3 = tftb.processing.cohen.smoothed_pseudo_wigner_ville(hilbert(x))
	#plt.figure()
	#plt.imshow(abs(wvd3), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
	#plt.xlabel('Time')
	#plt.ylabel('Frequency')
	#plt.colorbar()
	#plt.show()

