# -*- coding: utf-8 -*-
import pylab as py
import pickle
import numpy as np
import tftb
import matplotlib.pyplot as plt   
from SaveLoadData import *

def TFR(name):
	t,dt,Field,inc,dif,x,y,xmp,ymp,xm,ym = loadData(name)
	py.figure()
	py.plot(t,np.real(dif))    
	py.figure()
	py.scatter(x,y)  
	py.scatter(xm,ym)    
	py.show()
	tfr = tftb.processing.ShortTimeFourierTransform(np.real(dif),n_fbins = 4000,fwindow = np.hamming(2000))		#n_fbins
	tfr.run()
	tfr.plot(kind='cmap', show_tf=True)
	return tfr.tfr,tfr.ts,tfr.freqs

def quickStart():
	#TFR('trans2.pckl')
	#TFR('trans1.pckl')
	TFR('rot1.pckl')


if __name__ == '__main__':
	quickStart()
