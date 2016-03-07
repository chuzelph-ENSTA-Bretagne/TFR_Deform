import scipy.interpolate as inter
from scipy.signal import chirp, hilbert,find_peaks_cwt
import numpy as np
import pylab as plt
import tftb

import matplotlib.pyplot as plt
from skimage import data
try:
    from skimage import filters
except ImportError:
    from skimage import filter as filters
from skimage import exposure

def lawFrequency(tfr,t,f):
	tmp = np.zeros((len(f),len(t)))
	k = np.zeros(len(t))
	init = 0

	# creation du chemin
	for i in range(len(t)):
		if(i == 0 or init == 0):
			peakind = find_peaks_cwt(tfr[:,i], np.arange(1,10))
			if(len(peakind)!=0):
				print(peakind[0])
				tmp[peakind[0],i] = 1
				k[i] = peakind[0]
				init=1
			else:
				if(max(tfr[:,0])!=0):
					print('ici')
					a = 0
					var = 0
					for l in range(len(tfr[:,0])):
						if(tfr[l,i]>var):
							var = tfr[l,i]
							a = l
					tmp[a,i] = 1
					k[i]=a
					init=1
		else:
			Vartmp = 0
			Indice = 0
			for j in [-1,0,1]:
				if(Vartmp<tfr[k[i-1]+j,i]):
					Vartmp = tmp[k[i-1]+j,i]
					Indice = j
			k[i] = k[i-1] + Indice
			tmp[k[i],i] = 1
			#break
	#print(k)

	plt.figure()
	plt.imshow(tmp, origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
	plt.xlabel('Time')
	plt.ylabel('Frequency')
	plt.colorbar()
	plt.show()
	return k

#### Generation du signal
fm, am, iflaw = tftb.generators.misc.doppler(512, 200.0, 65.0, 10.0, 50.0)
sig = np.real(am * fm)
plt.figure("Doopler signal")
plt.subplot(211), plt.plot(sig)
plt.subplot(212), plt.plot(iflaw)






#### Obtention de la representation temps frequence
tfr = tftb.processing.ShortTimeFourierTransform(sig)
tfr.run()
tfr.plot(kind='cmap', show_tf=True)

wvd = tftb.processing.cohen.WignerVilleDistribution(hilbert(sig))
wvd.run()
wvd.plot(kind='cmap', show_tf=True) 




#### test des methodes de recuperation de la loi de frequence 




# methode : seuillage + seam carving

camera = np.abs(tfr.tfr)

# Filtrage d'otsu
val = filters.threshold_otsu(camera)

# representation de l'histogramme
hist, bins_center = exposure.histogram(camera)

# affichage des resultats
plt.figure(figsize=(9, 4))
plt.subplot(131)
plt.imshow(camera, cmap='jet', interpolation='nearest')
plt.axis('off')
plt.subplot(132)
plt.imshow(camera > val, cmap='jet', interpolation='nearest')
plt.colorbar()
plt.axis('off')
plt.subplot(133)
plt.plot(bins_center, hist, lw=2)
plt.axvline(val, color='k', ls='--')
plt.tight_layout()
plt.show()

VarSeuil = (camera > val)*tfr.tfr
plt.figure("VarSeuil signal")
plt.imshow(VarSeuil, cmap='jet', interpolation='nearest')
plt.colorbar()
plt.show()


lawFound = lawFrequency(np.abs(VarSeuil),tfr.ts,tfr.freqs)

#### spline sur la courbe trouvee

test = inter.UnivariateSpline(range(512), iflaw)
test.set_smoothing_factor(0.0005)
xnew = np.linspace(0, 512, num=512*10, endpoint=True)
plt.figure("spline pour la loi de frequence exact")
plt.plot (range(512), iflaw)
plt.plot (xnew, test(xnew),'+')
plt.show()


test2 = inter.UnivariateSpline(range(512), lawFound)
test2.set_smoothing_factor(0.0005)
xnew = np.linspace(0, 512, num=512*10, endpoint=True)
plt.figure("spline pour la loi de frequence trouvee")
plt.plot (range(512), lawFound)
plt.plot (xnew, test2(xnew),'+')
plt.show()
