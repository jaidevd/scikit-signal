import tf_generators as tfgen
import numpy as np
from matplotlib.pyplot import imshow, show, colorbar

N = 128
x = tfgen.fmlin(N, 0, 0.5)[0]
x = np.reshape(x,(len(x),1))

xrow, xcol = x.shape
t = np.arange(xrow)
trow, tcol = 1, len(t)


tfr = np.zeros((N, tcol), dtype=float)

for icol in range(1,tcol+1):
    ti = t[icol-1]
    taumax = min([ti-1, xrow-ti, round(N/2.0)-1])
    if icol == 1:
        tau = 0
    else:
        tau = np.arange(-taumax, taumax+1)
        
    indices = (N+tau)%N + 1
    
    tfr[indices-1, icol-1] = x[ti + tau -1, 0] * np.conj(x[ti-tau-1, xcol-1])
    tau = round(N/2.0)
    
    
    
    if (ti <= (xrow - tau)) or (ti >= (tau + 1)):
        try:
            tfr[tau, icol-1] = 0.5 * (x[ti+tau-1, 0]*np.conj(x[ti-tau-1, xcol-1]) + \
                                      x[ti-tau-1, 0]*np.conj(x[ti+tau-1, xcol-1]))
        except IndexError:
            pass

print tfr.shape


tfr = np.fft.fft(tfr, axis=0)
tfr = np.real(tfr)

imshow(tfr, aspect='auto')
colorbar()
show()