# -*- coding: utf-8 -*-
"""
Created on Thu Nov 05 13:31:06 2015

@author: Philippe
"""

import numpy as np
from scipy.signal import hilbert
    
def DWVD(s,fs,tau):
    """
    Discrete Wigner-Ville distribution based on paper:
    O'Toole, J., Mesbah, M., & Boashash, B. (2005). A discrete time and
    frequency Wigner-Ville distribution: properties and implementation.
    
    Input:
      s - signal 1D in callable array.
      
    Output:
      W - 2D numpy array of N x 2N size with frequencies on vertical axis.
    """
    
    # Obtenir le signal analytique
    # s = hilbert(s)    
    
    # Setting arrays
    N = s.size
    W = np.zeros((fs/2, N), dtype=np.complex64)
    
    piN = np.pi/float(N/tau)
    
    # For frequencies
    for k in range(fs/2):
        
        piNk = piN*k
    
        # For time shifts
        for n in range(N):
            
            # Limiting computational range
            l1 = max(0, n-N+1)
            l2 = min(n, N-1)
            
            m = np.arange(l1, l2)
            v = np.sum( s[m]*s[n-m]*np.exp(-2j*m*piNk) )
            
            W[k, n] = np.exp(1j*piNk*n)*v    
    
    return (np.real(W))    # np.flipud
    
