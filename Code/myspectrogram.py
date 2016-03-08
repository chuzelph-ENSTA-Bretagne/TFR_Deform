# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 15:14:25 2015

@author: Philippe
"""
import numpy as np
import scipy

# passage en python de la bibliotheque de Flandrin
# aide provenant egalement du site http://www.frank-zalkow.de/en/code-snippets/create-audio-spectrograms-with-python.html

def myspectrogram(x, wlen, h, nfft, fs):
    '''
    function: [stft, f, t] = stft(x, wlen, h, nfft, fs)
    x - signal in the time domain
    wlen - length of the hamming window
    h - hop size
    nfft - number of FFT points
    fs - sampling frequency, Hz
    stft - STFT matrix (only unique points, time across columns, freq across rows)
    Attention, nfft va améliorer la qualité en fréquence, h la qualité en temps
    '''
    # length of the signal
    xlen = np.size(x)
    
    # form a periodic hamming window
    win = scipy.hamming(wlen)
    
    # form the stft matrix
    rown = np.ceil((1+nfft)/2.0)             # calculate the total number of rows
    coln = 1+np.fix(float(xlen-wlen)/h)         # calculate the total number of columns
    stft = np.zeros([rown, coln],dtype=np.complex64)               # form the stft matrix
    
    # initialize the indexes
    indx = 0
    col = 1
    
    # perform STFT
    while (indx + wlen <= xlen):
        # windowing
        xw = x[indx:indx+wlen]*win;
        
        # FFT
        X = scipy.fft(xw, nfft);
        
        # update the stft matrix
        stft[:, col-1] = X[0:rown];
        
        # update the indexes
        indx = indx + h
        col = col + 1
    t = np.linspace(wlen/2,wlen/2+(coln+1)*h,coln)/fs
    f = np.linspace(0,fs/2,rown) #np.array(range(int(rown)))*fs/nfft
    return stft , t , f
    
    
