#!/usr/bin/env python

import numpy as np
from math import pi

def fmlin(N,fnormi=0.0,fnormf=0.5,t0=None):
    
    if not N:
        raise TypeError('Length of signal not provided.')
    
    if not t0:
        t0 = N/2
    
    if N <= 0:
        raise TypeError('Length of signal must be positive.')
    elif (abs(fnormi)>0.5) or (abs(fnormf>0.5)):
        raise TypeError('fnormi and fnormf must be between -0.5 and 0.5.')
    else:
        y = np.arange(N)
        y = fnormi*(y-t0) + ((fnormf - fnormi)/(2.0*(N-1))) * \
            ((y-1)**2 - (t0-1)**2)
        y = np.exp(1j*2.0*pi*y)
        y  = y/y[t0]
        
        iflaw = np.linspace(fnormi, fnormf, N)
    
    return y.T, iflaw