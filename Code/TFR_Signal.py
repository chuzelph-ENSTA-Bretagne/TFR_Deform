# -*- coding: utf-8 -*-
import pylab as py
import pickle
import numpy as np
import tftb
import matplotlib.pyplot as plt   
from SaveLoadData import *
from LawFrequency import *

def TFR(name,affichage1,affichage2,n_f,fSize):
	t,dt,Field,inc,dif,x,y,xmp,ymp,xm,ym = loadData(name)
	if(affichage1):
		py.figure()
		py.plot(t,np.real(dif))    
		py.figure()
		py.scatter(x,y)  
		py.scatter(xm,ym)    
		py.show()
	tfr = tftb.processing.ShortTimeFourierTransform(np.real(dif),n_fbins = n_f,fwindow = np.hamming(fSize))
	tfr.run()
	if(affichage2):
		tfr.plot(kind='cmap', show_tf=True)
	return tfr.tfr,tfr.ts,tfr.freqs

def TFR_WO_MEAN(name,affichage1,affichage2,n_f,fSize):
	t,dt,Field,inc,dif,x,y,xmp,ymp,xm,ym = loadData(name)
	if(affichage1):
		py.figure()
		py.plot(t,np.real(dif))    
		py.figure()
		py.scatter(x,y)  
		py.scatter(xm,ym)    
		py.show()
	tfr = tftb.processing.ShortTimeFourierTransform(np.real(dif)-np.mean(np.real(dif)),n_fbins = n_f,fwindow = np.hamming(fSize))
	tfr.run()
	if(affichage2):
		tfr.plot(kind='cmap', show_tf=True)
	return tfr.tfr,tfr.ts,tfr.freqs

def TFR_Sum(name1,name2,name3,affichage1,affichage2,n_f,fSize):
	t,dt,Field,inc,dif1,x,y,xmp,ymp,xm,ym = loadData(name1)
	t,dt,Field,inc,dif2,x,y,xmp,ymp,xm,ym = loadData(name2)
	t,dt,Field,inc,dif3,x,y,xmp,ymp,xm,ym = loadData(name3)
	if(affichage1):
		py.figure()
		py.plot(t,np.real(dif1)+np.real(dif2)+np.real(dif3))    
		py.figure()
		py.scatter(x,y)  
		py.scatter(xm,ym)    
		py.show()
	tfr = tftb.processing.ShortTimeFourierTransform(np.real(dif1)+np.real(dif2)+np.real(dif3),n_fbins = n_f,fwindow = np.hamming(fSize))		#n_fbins
	tfr.run()
	if(affichage2):
		tfr.plot(kind='cmap', show_tf=True)
	return tfr.tfr,tfr.ts,tfr.freqs

def TFR_Sum2(name1,name2,name3,affichage1,affichage2,n_f,fSize):
	t,dt,Field,inc,dif1,x,y,xmp,ymp,xm,ym = loadData(name1)
	t,dt,Field,inc,dif2,x,y,xmp,ymp,xm,ym = loadData(name2)
	t,dt,Field,inc,dif3,x,y,xmp,ymp,xm,ym = loadData(name3)
	if(affichage1):
		py.figure()
		py.plot(t,np.real(dif1)+np.real(dif2)+np.real(dif3))    
		py.figure()
		py.scatter(x,y)  
		py.scatter(xm,ym)    
		py.show()
	tfr1 = tftb.processing.ShortTimeFourierTransform(np.real(dif1),n_fbins = n_f,fwindow = np.hamming(fSize))	
	tfr2 = tftb.processing.ShortTimeFourierTransform(np.real(dif2),n_fbins = n_f,fwindow = np.hamming(fSize))	
	tfr3 = tftb.processing.ShortTimeFourierTransform(np.real(dif3),n_fbins = n_f,fwindow = np.hamming(fSize))
	tfr1.run()
	tfr2.run()
	tfr3.run()
	if(affichage2):
		tfr1.plot(kind='cmap', show_tf=True)
		tfr2.plot(kind='cmap', show_tf=True)
		tfr3.plot(kind='cmap', show_tf=True)
	return tfr1.tfr+tfr2.tfr+tfr3.tfr,tfr1.ts,tfr1.freqs

def quickStart():

	# RTF des signaux de base sans dilatation du rayon du cylindre			OK!
	#X,T,F = TFR('trans2.pckl',False,True,4000,2000)
	#X,T,F = TFR('trans1.pckl',False,True,4000,2000)
	#X,T,F = TFR('rot1.pckl',False,True,4000,2000)

	# Somme des TRF et RTF de la somme						OK!
	#X,T,F = TFR_Sum('trans1.pckl','trans2.pckl','rot1.pckl',False,False,4000,2000)
	#py.figure('RTF de la somme des trois signaux')
	#py.imshow(np.log10(abs(X)+1), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
	#py.xlabel('Time')
	#py.ylabel('Frequency')
	#py.colorbar()
	#py.show()

	#X,T,F = TFR_Sum2('trans1.pckl','trans2.pckl','rot1.pckl',False,False,4000,2000)
	#py.figure('Somme des RTF')
	#py.imshow(np.log10(abs(X)+1), origin='lower', cmap='jet', interpolation='nearest', aspect='auto')
	#py.xlabel('Time')
	#py.ylabel('Frequency')
	#py.colorbar()
	#py.show()

	# RTF du cylindre immobile et se contractant/dilatant				OK!
	#TFR_WO_MEAN('unanimate.pckl',False,True,4000,5000)


	# RTF des cylindres mobiles et se contractant/dilatant				OK!
	#TFR_WO_MEAN('trans1WD.pckl',False,True,4000,8000)



	# t,dt,Field,inc,dif,x,y,xmp,ymp,xm,ym = loadData('unanimate.pckl')			OK!
	# tfr = tftb.processing.ShortTimeFourierTransform(np.real(dif[0:len(dif)/100])-np.mean(np.real(dif[0:len(dif)/100])),n_fbins = len(dif)/100,fwindow = np.hamming(100))
	# tfr.run()
	# print(len(dif)/100)
	# py.figure()
	# py.plot(t[0:len(dif)/100],np.real(dif[0:len(dif)/100])) 
	# tfr.plot(kind='cmap', show_tf=True)

def verifSignal(n):
	for i in range(n+1):
		X,T,F = TFR("check"+str(i)+".pckl",False,True,4000,2000)


if __name__ == '__main__':
	#verifSignal(10)
	#quickStart()
	X,T,F = TFR('trans1.pckl',True,True,4000,2000)
	lawFound = seamCarving(X,T,F)
	splineFit(lawFound,T,0.0001)
