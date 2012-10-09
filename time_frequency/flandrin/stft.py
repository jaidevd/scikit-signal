#!/usr/bin/env python

import tf_generators as tfgen
import numpy as np
import matplotlib.pyplot as plt

f = tfgen.fmlin(128,0,0.5)[0]
window_lengths = [128, 64, 32, 16, 8, 4, 2, 1]

hamming_windows = []


for L in window_lengths:
    W = np.zeros((len(f),1))
    W[0:L] = np.hamming(L).reshape((L,1))
    stft = np.zeros((len(f)-L+1,len(f)), dtype=complex)
    for i in range(len(f)- L+1):
        indices = np.arange(i,i+L)
        x_stack = f*W
        stft[i,:] = x_stack.T[0]
        W = np.roll(W,1)
    hamming_windows.append(W)
    stft = np.real(stft)**2
    plt.subplot(2,4,window_lengths.index(L)+1)
    plt.imshow(stft[::-1], aspect = 'auto', cmap=plt.cm.jet)
    plt.title('L = '+str(L))


plt.show()

plt.figure()
for i in range(8):
    plt.subplot(2,4,i+1), plt.plot(hamming_windows[i])
plt.show()