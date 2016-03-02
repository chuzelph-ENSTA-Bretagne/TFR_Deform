# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 21:16:23 2015

@author: Philippe
MySpectrogram issu de Matlab
"""

import numpy as np
import scipy

def myspectrogram(x,nfft=2048,fs=1,window=scipy.hanning(512),noverlap=256,doplot=0,dbdown=100):

    '''
    %MYSPECTROGRAM Calculate spectrogram from signal.
    % B = MYSPECTROGRAM(A,NFFT,Fs,WINDOW,NOVERLAP) calculates the
    %     spectrogram for the signal in vector A.
    %
    % NFFT is the FFT size used for each frame of A.  It should be a
    % power of 2 for fastest computation of the spectrogram.
    %
    % Fs is the sampling frequency. Since all processing parameters are
    % in units of samples, Fs does not effect the spectrogram itself,
    % but it is used for axis scaling in the plot produced when
    % MYSPECTROGRAM is called with no output argument (see below).
    %
    % WINDOW is the length M window function applied, IN ZERO-PHASE
    % FORM, to each frame of A.  M cannot exceed NFFT.  For M<NFFT,
    % NFFT-M zeros are inserted in the FFT buffer (for interpolated
    % zero-phase processing).  The window should be supplied in CAUSAL
    % FORM.
    %
    % NOVERLAP is the number of samples the sections of A overlap, if
    % nonnegative.  If negative, -NOVERLAP is the "hop size", i.e., the
    % number of samples to advance successive windows.  (The overlap is
    % the window length minus the hop size.)  The hop size is called
    % NHOP below.  NOVERLAP must be less than M.
    %
    % If doplot is nonzero, or if there is no output argument, the
    % spectrogram is displayed.
    %
    % When the spectrogram is displayed, it is ``clipped'' dbdown dB
    % below its maximum magnitude.  The default clipping level is 100
    % dB down.
    %
    % Thus, MYSPECTROGRAM splits the signal into overlapping segments of
    % length M, windows each segment with the length M WINDOW vector, in
    % zero-phase form, and forms the columns of B with their zero-padded,
    % length NFFT discrete Fourier transforms.
    %
    % With no output argument B, MYSPECTROGRAM plots the dB magnitude of
    % the spectrogram in the current figure, using
    % IMAGESC(T,F,20*log10(ABS(B))), AXIS XY, COLORMAP(JET) so the low
    % frequency content of the first portion of the signal is displayed
    % in the lower left corner of the axes.
    %
    % Each column of B contains an estimate of the short-term,
    % time-localized frequency content of the signal A.  Time increases
    % linearly across the columns of B, from left to right.  Frequency
    % increases linearly down the rows, starting at 0.
    %
    % If A is a length NX complex signal, B is returned as a complex
    % matrix with NFFT rows and
    %      k = floor((NX-NOVERLAP)/(length(WINDOW)-NOVERLAP))
    %        = floor((NX-NOVERLAP)/NHOP)
    % columns.  When A is real, only the NFFT/2+1 rows are needed when
    % NFFT even, and the first (NFFT+1)/2 rows are sufficient for
    % inversion when NFFT is odd.
    %
    % See also: Matlab's SPECTROGRAM and Octave's STFT function.
    '''

    M = len(window);
    if (M<2):
        print('myspectrogram: Expect complete window, not just its length')
        
    if (len(x)<M):
        #zero-pad to fill a window:
        x = np.concatenate((x , np.zeros(M-len(x))))
    Modd = np.mod(M,2) # 0 if M even, 1 if odd
    Mo2 = (M-Modd)/2.0
    
    if noverlap<0:
      nhop = - noverlap
      noverlap = M-nhop
    else:
      nhop = M-noverlap
    
    nx = len(x)
    #nframes = 1+floor((nx-noverlap)/nhop);
    nframes = int(1+np.ceil(float(nx)/float(nhop)));
    
    X = np.zeros([nfft,nframes]) # allocate output spectrogram
    
    zp = np.zeros(nfft-M)  #zero-padding for each FFT
    xframe = np.zeros(M);
    xoff = 0 - Mo2 # input time offset = half a frame
    for m in range(nframes):
    #  M,Mo2,xoff,nhop
        if (xoff<0):
            xframe[0:xoff+M] = x[0:xoff+M] # partial input data frame
        else:
            if (xoff+M > nx):
                xframe = np.concatenate((x[xoff:nx],np.zeros(xoff+M-nx)))
            else:
                xframe = x[xoff:xoff+M] # input data frame
        xw = window*xframe # Apply window
        xwzp = np.concatenate((xw[Mo2:M],zp,xw[0:Mo2]))
        X[:,m] = scipy.fft(xwzp)
        xoff = xoff + nhop # advance input offset by hop size
    t = np.linspace(0,nframes-1)*nhop/fs;
    f = 0.001*np.linspace(0,nfft-1)*fs/nfft;
    return X,t,f
    
#    if (nargout==0) | doplot
#      t = (0:nframes-1)*nhop/fs;
#      f = 0.001*(0:nfft-1)*fs/nfft;
#      Xdb = 20*log10(abs(X));
#      Xmax = max(max(Xdb));
#      % Clip lower limit to -dbdown dB so nulls don't dominate:
#      clipvals = [Xmax-dbdown,Xmax];
#      imagesc(t,f,Xdb,clipvals);
#      % grid;
#      axis('xy');
#      colormap(jet);
#      xlabel('Time (sec)');
#      ylabel('Freq (kHz)');
#    end