# -*- coding: utf-8 -*-

import numpy as np
import math as m
import scipy
import pylab as py
from GenerationSignalV2 import *
import pickle
from testSave import *


class ObjectDiff():
    
    '''
    x,y                 : position de l'objet
    theta,rayon             : coordonnées polaire de l'objet, calcule a chaque instant
    alpha               : parametre du probleme, angle de l'onde incidente
    xm, ym              : position de l'observateur, supposé constant dans toute la simulation
    thetaDif , rDiff    : distance du centre de l'objet à l'observateur et son angle
    '''
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
        
    def runOneStep(self):
        self.updatePos()
        self.phase = self.x*self.kx+self.y*self.ky
        self.Field_incTMP   = np.exp(1j*(self.xm*self.kx+self.ym*self.ky))
        self.Field_scatTMP  = Field_scat(self.k,self.r,self.rp,self.alpha,self.thetap,self.N)*np.exp(1j*self.phase)
        self.Field  = self.Field_incTMP + self.Field_scatTMP
        self.iter = self.iter +1
        
    def updatePos(self):
		xn,yn = self.Deplacement(self.x,self.y)
		self.x = xn 
		self.y = yn
		self.xmp = self.xm-self.x
		self.ymp = self.ym-self.y
		self.rp = np.sqrt((self.xmp)**2+(self.ymp)**2)
		self.thetap = np.arctan2((self.ymp),(self.xmp))
		self.r = self.radius0 + self.r0var*np.sin(2*m.pi*self.freq*self.t[self.iter])
    
    def Deplacement(self,x,y):
        X = np.array([x,y])
        c = np.cos(self.Vtheta*self.dt)
        s = np.sin(self.Vtheta*self.dt)
        R = np.array([[c,-s],[s,c]])
        X = np.dot(R,X) + self.dt*np.array([self.Vx , self.Vy])
        return X[0],X[1]
        
    def affichage(self):
        print("x :",self.x)
        print("y :",self.y)
        print("xm :",self.xm)
        print("ym :",self.ym)
        print("r :",self.r)
        print("rDiff :",self.rDiff)
        print("thetaDif :",self.thetaDif)
        print("N :",self.N)
        
        
        return self.x,self.y,self.xm,self.ym
        
        
if __name__ == '__main__': 

	k=5
	radius0=1.
	alpha=0.
	x=0
	y=10
	xm = 10.
	ym = 0.


	freq=250.
	tau = 4
	nbpts = 10*tau*freq
	t = np.linspace(0,tau,nbpts)
	dt = tau/nbpts


	Obj1 = ObjectDiff(t,dt,x,y,alpha,xm,ym,k,radius0=1,r0var=0.0,freq = 250,Vx = 10,Vy = 0, Vtheta = 0*2*m.pi/360)

	List = []
	x = []
	y = []
	phase = []
	inc = []
	dif = []
	rp = []
	thetap = []
	xmp = []
	ymp = []
	for i in range(len(t)):
		Obj1.runOneStep()
		List.append(Obj1.Field)
		x.append(Obj1.x)
		y.append(Obj1.y)
		phase.append(Obj1.phase)
		inc.append(Obj1.Field_incTMP)
		dif.append(Obj1.Field_scatTMP)
		rp.append(Obj1.rp)
		thetap.append(Obj1.thetap)
		xmp.append(Obj1.xmp)
		ymp.append(Obj1.ymp)
        
	List2 = np.array(List)
	List3 = np.array(inc)
	List4 = np.array(dif)

	py.figure('signal obtenu')
	py.plot(t,np.real(List2))  
	py.plot(t,np.real(List3))  
	py.plot(t,np.real(List4))    
	py.xlabel('time')
	py.ylabel('Real value of the signal')

	py.figure('position')
	py.scatter(x,y)
	py.scatter(Obj1.xm,Obj1.ym)
	py.xlabel('x')
	py.ylabel('y')

	py.figure('phase en cours du temps')
	py.plot(t,np.array(phase))
	py.xlabel('time')
	py.ylabel('difference de marche')

	py.figure('distance a emmeteur')
	py.plot(t,np.array(rp))
	py.xlabel('time')
	py.ylabel('distance a emmeteur')

	py.figure('angle avec emmeteur')
	py.plot(t,np.array(thetap)*180/(m.pi))
	py.xlabel('time')
	py.ylabel('angle avec emmeteur')

	#py.figure('xmp')
	#py.plot(t,np.array(xmp))
	#py.xlabel('time')
	#py.ylabel('xmp')

	#py.figure('ymp')
	#py.plot(t,np.array(ymp))
	#py.xlabel('time')
	#py.ylabel('ymp')
	py.show()

	saveData(t,np.real(List2),x,y,'Tran2.pckl')

    
    
    
