#!/usr/bin/env python

import numpy as np

def nextpow2(X):
    N = 0
    while 2**N < abs(X):
        N += 1
    return N

def tfrstft(x,t=None,N=None,h=None,trace=0):
    
    if t is None:
        t = np.arange(1,len(x)+1)
    
    if not N:
        N = len(x)
    
    if not h:
        h = np.hamming(N/4)
    
    xrow, xcol = x.shape
    
    hlength = np.floor(N/4)
    hlength = hlength + 1 - hlength%2
    
    try:
        trow, tcol = t.shape
    except ValueError:
        trow = 1
        tcol = len(t)
    
    if xcol!=1:
        raise IndexError("X must have only one column.")
    elif trow!=1:
        raise IndexError("T must have only one row.")
    
    
    try:
        hrow = h.shape[0]
    except AttributeError:
        hrow = 1
    Lh = (hrow-1)/2

    h = h/np.linalg.norm(h)
        
    tfr = np.zeros((N,tcol))
    
    
    if trace:
        print "Short Time Fourier Transform\n"
    
    for icol in range(1,tcol+1):
        ti = t[icol-1]
        taumin = min([round(N/2)-1, Lh, ti-1])
        taumax = min([round(N/2)-1, Lh, xrow-ti])
        tau = np.arange(-taumin, taumax+1)
        indices = (N+tau)%N +1
        tfr[indices-1, icol-1] = x[ti+tau-1, 0]
    tfr =  np.fft.fft(tfr, axis=0)
    
    return tfr

def tfrqview(tfr, sig=np.array([]), t=None, method='type1', **kwargs):
    
    [tfrrow, tfrcol] = tfr.shape
    N = tfrcol
    
    if not t:
        t = np.arange(N)
    
    trow, tcol = 1, len(t)
    Nsig, Ncol = sig.shape
    
    if method == 'type2':
        Nf2 = tfrrow;
    else:
        Nf2 = tfrrow/2
    
    freq = 0.5*np.arange(Nf2-2)/Nf2
    
    
    # Visualization options
    ###########################################################################
    threshold = 0.5
    linlogtfr = 0
    linlogspec = 1
    sigenveloppe = 0
    
    levelnumb = 64
    colmap = 1
    
    display = 2
    
    isgridsig = 0
    isgridspec = 0
    isgridtfr = 0
    
    issig = 0
    isspec = 0
    iscolorbar = 0
    
    fs = 1.0
    fmin = 0.0
    fmax = 0.5*fs
    
    
    # Test of Analyticity
    ###########################################################################
    if sig.shape == (0,):
        
        for k in range(Ncol):
            Lt_fog = np.max(t) - np.min(t) +1
            Nb_tranches_fog = np.floor(Lt_fog/tfrrow)
            spec[:,k-1] = np.zeros((tfrrow,1), dtype=float)
            
            for Num_tranche_fog in range(Nb_tranches_fog):
                fftarg = sig[
                    min(t)+tfrrow*Num_tranche_fog + np.arange(tfrrow) - 1,
                    k-1
                ]
                spec[:,k-1] = spec[:,k-1] + abs(np.fft.fft(fftarg,axis=0))**2
            
            if (Lt_fog>Nb_tranches_fog*tfrrow):
                #fftarg = sig(min(t)+tfrrow*Nb_tranches_fog:max(t),k),tfrrow
                fftarg = sig[
                    np.amin(t)+np.arange(tfrrow*Nb_tranches_fog, np.amax(t)+1)-1,
                    k-1
                ]
                spectre_fog = np.fft.fft(fftarg, n=tfrrow, axis=0)
                spectre_fog = spectre_fog.reshape((spectre_fog.size,1))
                spec[:,k-1] = spec[:,k-1] +abs(spectre_fog)**2
            
            spec1 = sum(spec[0:tfrrow/2, k-1])
            spec2 = sum(spec[tfrrow/2:tfrrow, k-1])
            
            if spec2 > spec1/10:
                print "Caution: the signal is not analytic!\n"
    
    
    
    # Test of Reality
    ###########################################################################
    if (Ncol==2) and np.isreal(tfr).all():
        print "Cross distribution. Result is complex, displaying real part.\n"
        tfr = np.real(tfr)
    
    ChoiceDisplay = 1
    ChoiceLayout      =  2
    ChoiceSampling    =  3
    ChoiceFreqBounds  =  4
    ChoiceThreshold   =  5
    ChoiceLinlog      =  6
    ChoiceRedraw      =  7
    ChoiceNewFigure   =  8
    ChoiceSaveResults =  9
    ChoiceSaveOptions = 10
    ChoicePrint       = 11
    ChoiceClose       = 12
    
    CallTfrView = 1
    RefreshFigure = 1
    choice = ChoiceSampling
    while choice != ChoiceClose:
        if RefreshFigure and CallTfrView:
            linlog = linlogtfr + 2*linogspec + 4*sigenveloppe
            isgrid = isgridsig*2 + iscolorbar*4 + 1
            param = [
                display, linlog, threshold, levelnumb, Nf2, layout, fs, isgrid,
                fmin, fmax
            ]
            

def tfrsp(x, t=None, N=None, h=None, trace=0):
    
    if not N:
        N = len(x)
    
    if not t:
        t = np.arange(N)
    
    if not h:
        h = np.hamming(N/4).reshape((N/4,1))
    
    xrow, xcol = x.shape
    hlength = np.floor(N/4)
    hlength = hlength + 1 - hlength%2
    trow, tcol = 1, len(t)
    hrow, hcol = h.shape
    Lh = (hrow-1)/2
    
    tfr = np.zeros((N, tcol), dtype=float)
    
    for i in range(tcol):
        ti = t[i]
        left = min([round(N/2)-1, Lh, ti-1])
        right = min([round(N/2)-1, Lh, xrow-ti])
        tau = np.arange(-left, right+1)
        indices = (N+tau)%N + 1
        rhs = x[ti+tau-1]*np.conj(h[Lh+tau])/np.linalg.norm(h[Lh+tau], ord=2)
        tfr[indices-1, i] = rhs.T[0]
    
    tfr = abs(np.fft.fft(tfr, axis=0))**2
    
    return tfr

def tfrwv(x, t=None, N=None, trace=0):
    
    if not x:
        raise TypeError('Atleast one parameter required.')
    
    if x.ndim == 1:
        xrow, xcol = 1, x.shape[0]
    else:
        xrow, xcol = x.shape
    
    if not t:
        t = np.arange(1,xrow+1)
    
    if not N:
        N = xrow
    
    else:
        if N < 0:
            raise TypeError('N must be greater than zero.')
    
    trow, tcol = t.shape
    
    if (xcol==0) or (xcol>2):
        raise TypeError('X must have one or two columns')
    elif trow!=1:
        raise TypeError('T must have only one row')
    elif (2**nextpow2(N)!=N):
        print 'For faster computation N should be a power of two.\n'
    
    tfr = np.zeros((N, tcol), dtype=float)
    
    if trace:
        print 'Wigner Ville distribution'
    
    for icol in range(tcol):
        ti = t[icol]
        taumax = min()