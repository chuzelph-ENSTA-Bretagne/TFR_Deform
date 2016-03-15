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




#### fonction qui prend comme parametre d'entree les sorties de la TFR et renvoie en sortie la premiere loi de frequence detectee
#### ATTENTION, il faut que tfr soit en valeur absolue pour que l'algorithme suivant converge
def lawFrequency(tfr,t,f):

	# Image qui sert a voir le resultat de la fonction
	tmp = np.zeros((len(f),len(t)))

	# Variable qui va stocker les valeurs trouvees
	k = np.zeros(len(t))

	# variable temporaire qui permet de savoir si une loi a ete trouve au debut de la fonction
	init = 0

	# creation du chemin en realisant une boucle sur le temps
	for i in range(len(t)):

		# phase d'initialisation, detecte le premier pic visible ou le premier maximum non nul le cas echeant
		if(i == 0 or init == 0):

			# find_peaks_cwt de la bibliotheque scipy.signal, assez pratique pour trouver tous les maximums d'une fonction. Peut etre utiliser pour trouver plusieurs maximum et suivre plusieurs loi en meme temps.
			peakind = find_peaks_cwt(tfr[:,i], np.arange(1,10))

			# petit test pour voir si un pic a bien ete trouve. Si oui, on stocke le coefficient trouve et on passe dans le cas else
			if(len(peakind)!=0):
				tmp[peakind[0],i] = 1
				k[i] = peakind[0]
				init=1



			# test si un maximum a pu etre trouve dans le cas ou la fonction find_peaks_cwt n'a pas ete concluante
			else:
				if(max(tfr[:,i])!=0):
					a = 0
					var = 0
					# boucle qui sert a trouve la position du max... Cela n'est pas a priori une fonction par defaut de python...
					for l in range(len(tfr[:,i])):
						if(tfr[l,i]>var):
							var = tfr[l,i]
							a = l
					tmp[a,i] = 1
					k[i]=a
					init=1

		# Une premiere loi a ete trouve et on construit pas a pas le chemin suivant l'algorithme du seams carving
		else:
			Vartmp = 0
			Indice = 0
			for j in [-1,0,1]:
				if(Vartmp<tfr[k[i-1]+j,i]):
					Vartmp = tfr[k[i-1]+j,i]
					Indice = j
			k[i] = k[i-1] + Indice
			tmp[k[i],i] = 1


	# Affiche la figure qui represent le chemin retenu
	plt.figure()
	plt.imshow(tmp, origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
	plt.xlabel('Time')
	plt.ylabel('Frequency')
	plt.colorbar()

	# renvoie la loi frequence en fonction du temps
	return k


#### methode : seuillage + seam carving
def seamCarving(tfr,t,f):

	# calcul la valeur absolue de l'image
	Image = np.abs(tfr)
	# Filtrage d'otsu
	val = filters.threshold_otsu(Image)

	# representation de l'histogramme
	hist, bins_center = exposure.histogram(Image)

	# affichage des resultats
	plt.figure(figsize=(9, 4))
	plt.subplot(131)
	plt.imshow(Image, cmap='jet', interpolation='nearest')
	plt.axis('off')
	plt.subplot(132)
	plt.imshow(Image > val, cmap='jet', interpolation='nearest')
	plt.colorbar()
	plt.axis('off')
	plt.subplot(133)
	plt.plot(bins_center, hist, lw=2)
	plt.axvline(val, color='k', ls='--')
	plt.tight_layout()

	# Creation du masque qui garde que les points significatifs
	VarSeuil = (Image > val)*tfr

	# affichage du mask
	plt.figure("VarSeuil signal")
	plt.imshow(np.flipud(VarSeuil), cmap='jet', interpolation='nearest')
	plt.colorbar()

	# fait appel a la fonction lawFrequency qui implemente le seam carving
	lawFound = lawFrequency(np.abs(VarSeuil),t,f)
	plt.show()
	return lawFound


#### Methode qui realise un fit spline sur la courbe retenue. ATTENTION signal 1D en entree
def splineFit(lawFound,t,smooth):

	# On doit en premier lieu rogner le signal afin d'enlever des zeros au debut de la loi si l'algorithme du seams carving n'a pas correctement fonctionne au debut du programme.
	tmp = 0
	while(lawFound[0]==0):
		tmp=tmp+1
		lawFound=lawFound[1:len(lawFound)]
	# Interpolation
	test = inter.UnivariateSpline(range(len(t)-tmp), lawFound)

	# lissage
	test.set_smoothing_factor(smooth)

	# definition d'un interval de temps avec plus de point
	Tnew = np.linspace(t[tmp], t[len(t)-1], num=(len(t)-tmp), endpoint=True)

	# Affichage du resultat
	plt.figure("spline pour la loi de frequence trouvee")
	plt.plot (t[tmp:len(t)], lawFound)
	plt.plot (Tnew, test(Tnew),'+')
	plt.show()



# test unitaire
if __name__ == '__main__':

	
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

	#### Appel des fonctions
	lawFound = seamCarving(tfr.tfr,tfr.ts,tfr.freqs)
	splineFit(lawFound,tfr.ts,0.0001)
