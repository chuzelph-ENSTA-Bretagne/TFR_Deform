# -*- coding: utf-8 -*-

import numpy as np
import math as m
import scipy
import pylab as py
from GenerationSignal import *
import pickle
from SaveLoadData import *

#### ATTENTION!!!!!!!!!!!
#### Si l'objet et l'observateur se trouve sur la meme position, L'ENSEMBLE DU CALCUL SERA INEXPLOITABLE!!!

class ObjectDiff():

	'''
	x,y                 : position de l'objet
	xm, ym              : position de l'observateur, supposé constant dans toute la simulation
	xmp, ymp            : position de l'observateur dans le repere du cylindre
	rp, thetap          : position de l'observateur dans le repere du cylindre en coordonnees polaire
	alpha               : parametre du probleme, angle de l'onde incidente
	radius0,r0var		: rayon initial et taux de variation du rayon du cylindre
	freq				: fréquence de dillatation/contraction du cylindre
	Vx,Vy				: vitesse de l'objet dans l'axe Ox et Oy
	Vtheta				: vitesse de rotation de l'objet selon l'axe Oz
	k					: paramètre de l'onde
	t,dt				: interval de temps de la simulation et son pas, a donner ABSOLUMENT au debut de la simulation
	'''


	#### INITIALISATION DE L'OBJET ####
	def __init__(self,t,dt,x,y,alpha,xm,ym,k,radius0=1,r0var=0.1,freq = 0,Vx = 0,Vy = 0, Vtheta = 0):

		# initialisation de l'objet
		self.t = t
		self.dt = dt
		self.iter = 0

		# parametre geometrique du probleme
		self.x = x
		self.y = y
		self.xm = xm
		self.ym = ym
		self.alpha = alpha
		self.xmp = xm-x
		self.ymp = ym-y
		self.rp = np.sqrt((self.xmp)**2+(self.ymp)**2)
		self.thetap = np.arctan2((self.ymp),(self.xmp))

		# Variation du rayon
		self.radius0 = radius0
		self.r = radius0
		self.r0var = r0var
		self.freq = freq

		# paramètre pour le déplcament de l'objet
		self.Vx = Vx				# m.s-1
		self.Vy = Vy				# m.s-1
		self.Vtheta = Vtheta		# rad.s-1
		self.phase = 0

		# parametres pour la generation de l'onde
		self.Field = 0
		self.Field_incTMP = 0
		self.Field_scatTMP = 0
		self.k = k
		self.kx = -np.cos(self.alpha)*k
		self.ky = np.sin(self.alpha)*k
		self.N = int(np.ceil(k*self.radius0*(self.radius0+self.r0var) + 10))
		

	#### Realise un tour de boucle ####
	def runOneStep(self):
		self.updatePos()
		self.phase = self.x*self.kx+self.y*self.ky
		self.Field_incTMP   = np.exp(1j*(self.xm*self.kx+self.ym*self.ky))
		self.Field_scatTMP  = Field_scat(self.k,self.r,self.rp,self.alpha,self.thetap,self.N)*np.exp(1j*self.phase)
		self.Field  = self.Field_incTMP + self.Field_scatTMP
		self.iter = self.iter +1
		
	#### Met a jour la position de l'objet ####
	def updatePos(self):
		xn,yn = self.Deplacement(self.x,self.y)
		self.x = xn 
		self.y = yn
		self.xmp = self.xm-self.x
		self.ymp = self.ym-self.y
		self.rp = np.sqrt((self.xmp)**2+(self.ymp)**2)
		self.thetap = np.arctan2((self.ymp),(self.xmp))
		self.r = self.radius0 + self.r0var*np.sin(2*m.pi*self.freq*self.t[self.iter])

	#### Calcul la nouvelle position de l'objet ####
	def Deplacement(self,x,y):
		X = np.array([x,y])
		c = np.cos(self.Vtheta*self.dt)
		s = np.sin(self.Vtheta*self.dt)
		R = np.array([[c,-s],[s,c]])
		X = np.dot(R,X) + self.dt*np.array([self.Vx , self.Vy])
		return X[0],X[1]
		
	#### affichage de la position de l'objet, test unitaire ####
	def affichage(self):
		print("x :",self.x)
		print("y :",self.y)
		print("xm :",self.xm)
		print("ym :",self.ym)
		print("xmp :",self.xmp)
		print("ymp :",self.ymp)
		print("N :",self.N)
		return None


	def run(self,affichage,save,nom):
		Field = []
		x = []
		y = []
		phase = []
		inc = []
		dif = []
		rp = []
		thetap = []
		xmp = []
		ymp = []
		for i in range(len(self.t)):
			self.runOneStep()
			Field.append(self.Field)
			x.append(self.x)
			y.append(self.y)
			phase.append(self.phase)
			inc.append(self.Field_incTMP)
			dif.append(self.Field_scatTMP)
			rp.append(self.rp)
			thetap.append(self.thetap)
			xmp.append(self.xmp)
			ymp.append(self.ymp)

		if(affichage):
			py.figure('signal obtenu')
			py.plot(self.t,np.real(np.array(Field)))  
			py.plot(self.t,np.real(np.array(inc)))  
			py.plot(self.t,np.real(np.array(dif)))    
			py.xlabel('time')
			py.ylabel('Real value of the signal')

			py.figure('position')
			py.scatter(x,y)
			py.scatter(self.xm,self.ym)
			py.xlabel('x')
			py.ylabel('y')

			py.figure('difference de marche au cours du temps')
			py.plot(self.t,np.array(phase))
			py.xlabel('time')
			py.ylabel('difference de marche')

			py.figure('distance a emmeteur')
			py.plot(self.t,np.array(rp))
			py.xlabel('time')
			py.ylabel('distance a emmeteur')

			py.figure('angle avec emmeteur')
			py.plot(self.t,np.array(thetap)*180/(m.pi))
			py.xlabel('time')
			py.ylabel('angle avec emmeteur')

			py.figure('xmp')
			py.plot(self.t,np.array(xmp))
			py.xlabel('time')
			py.ylabel('xmp')

			py.figure('ymp')
			py.plot(self.t,np.array(ymp))
			py.xlabel('time')
			py.ylabel('ymp')
			py.show()

		if(save):
			saveData(self.t,self.dt,np.array(Field),np.array(inc),np.array(dif),np.array(x),np.array(y),np.array(xmp),np.array(ymp),self.xm,self.ym,nom)

#### Petite metthode qui permet de generer des signaux simples lorsque vous charger le git.
#### Peut servir de test unitaire pour voir comment utiliser la classe generee
def quickStart():
	# parametres du probleme
	k=100
	radius0=1.
	alpha=0.
	x=0
	y=10
	xm = 5.
	ym = 0.
	freq=500.							# frequence de dillatation/contraction du cylindre s-1
	tau = 2								# duree de la simulation en seconde
	nbpts = 20*tau*max(freq,k)					# nombre de point souhaite. PENSER A METTRE BEAUCOUP DE POINT!!!!!!!
	t = np.linspace(0,tau,nbpts)		# Interval de temps d'etude. A donner en debut de simulation
	dt = tau/nbpts						# pas de temps, a calculer au debut du programme


	# Creation de trois objets, deux en translation, un en rotation.
	print("trans1.pckl...")
	Obj1 = ObjectDiff(t,dt,x,y,alpha,xm,ym,k,radius0=1,r0var=0.0,freq = 100,Vx = 5,Vy = 0, Vtheta = 0*2*m.pi/360)
	Obj1.run(False,True,"trans1.pckl")
	print("Done")
	print("trans2.pckl...")
	Obj2 = ObjectDiff(t,dt,10,y,alpha,xm,ym,k,radius0=1,r0var=0.0,freq = 100,Vx = -5,Vy = 0, Vtheta = 0*2*m.pi/360)
	Obj2.run(False,True,"trans2.pckl")
	print("Done")
	print("rot1.pckl...")
	Obj3 = ObjectDiff(t,dt,3,0,alpha,xm,ym,k,radius0=1,r0var=0.0,freq = 100,Vx = 0,Vy = 0, Vtheta = 90*2*m.pi/360)
	Obj3.run(False,True,"rot1.pckl")
	print("Done")

	# Creation d'un objet immobile, mais qui se dilate/contracte.
	print("unanimate.pckl...")
	Obj3 = ObjectDiff(t,dt,3,0,alpha,xm,ym,k,radius0=1,r0var=0.1,freq = 100,Vx = 0,Vy = 0, Vtheta = 0*2*m.pi/360)
	Obj3.run(False,True,"unanimate.pckl")
	print("Done")

	# Creation de trois objets, deux en translation, un en rotation, tous avec une dilatation / contraction du rayon.
	print("trans1WD.pckl...")
	Obj3 = ObjectDiff(t,dt,x,y,alpha,xm,ym,k,radius0=1,r0var=0.1,freq = 100,Vx = 5,Vy = 0, Vtheta = 0*2*m.pi/360)
	Obj3.run(False,True,"trans1WD.pckl")
	print("Done")
	print("trans2WD.pckl...")
	Obj3 = ObjectDiff(t,dt,10,y,alpha,xm,ym,k,radius0=1,r0var=0.1,freq = 100,Vx = 0,Vy = -5, Vtheta = 0*2*m.pi/360)
	Obj3.run(False,True,"trans2WD.pckl")
	print("Done")
	print("rot1WD.pckl...")
	Obj3 = ObjectDiff(t,dt,3,0,alpha,xm,ym,k,radius0=1,r0var=0.1,freq = 100,Vx = 0,Vy = 0, Vtheta = 90*2*m.pi/360)
	Obj3.run(False,True,"rot1WD.pckl")
	print("Done")


if __name__ == '__main__': 
	quickStart()
	print("Script ObjectDiff.py termine")
    
