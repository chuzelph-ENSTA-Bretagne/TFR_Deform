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
    theta,r             : coordonnées polaire de l'objet, calcule a chaque instant
    alpha               : parametre du probleme, angle de l'onde incidente
    xm, ym              : position de l'observateur, supposé constant dans toute la simulation
    thetaDif , rDiff    : distance du centre de l'objet à l'observateur et son angle
    '''
    def __init__(self,t,dt,x,y,alpha,xm,ym,k,radius0=1,r0var=0.1,freq = 0,Vt = 0,Thetat = 0, Vr = 0):
        
        # initialisation de l'objet
        self.t = t
        self.dt = dt
        self.iter = 0
        
        # parametre geometrique du probleme
        self.x = x
        self.y = y
        self.xm = xm
        self.ym = ym
        self.rm = np.sqrt((xm)**2+(ym)**2)
        self.Thetam = np.arctan2((ym),(xm))
        self.alpha = alpha
        self.rDiff = np.sqrt((xm-x)**2+(ym-y)**2)
        self.thetaDif = np.arctan2((ym-y),(xm-x))
        
        # Variation du rayon
        self.radius0 = radius0
        self.r = radius0
        self.r0var = r0var
        self.freq = freq
        
        # paramètre pour le déplcament de l'objet
        self.Vt = Vt
        self.Thetat = Thetat
        self.Vr = Vr
        self.phase = 0
        
        # parametres pour la generation de l'onde
        self.Field = 0
        self.k = k
        self.N = int(np.ceil(k*self.radius0*(self.radius0+self.r0var) + 10))
        
    def runOneStep(self):
        self.updatePos()
        rayon = np.sqrt((self.x)**2+(self.y)**2)
        theta = np.arctan2((self.y),(self.x))
        self.phase = rayon*np.cos(self.alpha-theta)
        Field_incTMP   = np.exp(1j*self.k*self.rm*np.cos(self.Thetam-self.alpha))*np.exp(1j*k*self.phase)
        Field_scatTMP  = Field_scat(self.k,self.r,self.rDiff,self.alpha,self.thetaDif,self.N)
        self.Field  = Field_incTMP + Field_scatTMP
        self.iter = self.iter +1
        
    def updatePos(self):
        xn,yn = self.Deplacement(self.x,self.y)
        self.x = xn 
        self.y = yn
        self.rDiff = np.sqrt((xm-xn)**2+(ym-yn)**2)
        self.thetaDif = np.arctan2((ym-yn),(xm-xn))
        self.r = self.radius0 + self.r0var*np.sin(2*m.pi*self.freq*self.t[self.iter])
    
    def Deplacement(self,x,y):
        X = np.array([x,y])
        c = np.cos(self.Vr*self.dt)
        s = np.sin(self.Vr*self.dt)
        R = np.array([[c,-s],[s,c]])
        X = np.dot(R,X) + self.Vt*self.dt*np.array([np.cos(self.Thetat) , np.sin(self.Thetat)])
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
    r=10.
    theta=0.
    radius0=1.
    alpha=0.
    x=-15
    y=3
    xm = r*np.cos(theta)
    ym = r*np.sin(theta)
    
    
    freq=250.
    tau = 10
    nbpts = 10*tau*freq
    t = np.linspace(0,tau,nbpts)
    dt = tau/nbpts
    
    
    Obj1 = ObjectDiff(t,dt,x,y,alpha,xm,ym,k,radius0=1,r0var=0.0,freq = 250,Vt = 5,Thetat = 0, Vr = 0*2*m.pi/360)
    
    List = []
    x = []
    y = []
    phase = []
    for i in range(len(t)):
        Obj1.runOneStep()
        List.append(Obj1.Field)
        x.append(Obj1.x)
        y.append(Obj1.y)
        phase.append(Obj1.phase)
        
    List2 = np.array(List)
    py.figure()
    py.plot(t,np.real(List2))    
    py.figure()
    py.scatter(x,y)
    py.scatter(Obj1.xm,Obj1.ym)
    py.figure()
    py.plot(t,np.array(phase))
    py.show()
    
    saveData(t,np.real(List2),x,y,'Tran2.pckl')
    
    
    
    