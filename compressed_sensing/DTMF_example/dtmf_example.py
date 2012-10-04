#!/usr/bin/env python

import pyglet
import DTMF
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model.sparse import Lasso
from sklearn.linear_model import Ridge

FR = 44000
scale = 32767
keys=   '1','2','3','A',\
	'4','5','6','B',\
	'7','8','9','C',\
	'*','0','#','D'

F1 = [697,770,852,941]
F2 = [1209, 1336, 1477, 1633]


class ArraySource(pyglet.media.StaticMemorySource):
    '''A source that has been created from a numpy array.'''

    def __init__(self, rate, data):
        '''Construct an `ArraySource` for the data in `data`.

        :Parameters:
            `rate` : `int`
                The sampling rate in Hertz.
            `data` : `numpy.array`
                A c-contiguous numpy array of dimension (`num_samples`,
                `num_channels`).
        '''

        if data.ndim not in (1, 2):
            raise ValueError("The data array must be one- or two-dimensional.""")

        if not data.flags.c_contiguous:
            raise ValueError("The data array must be c-contiguous.""")

        num_channels = data.shape[1]
        if num_channels not in [1, 2]:
            raise ValueError("Only mono and stereo audio are supported.""")

        num_bits = data.dtype.itemsize
        if num_bits not in [8, 16]:
            raise ValueError("Only 8 and 16 bit audio are supported.""")

        audio_format = pyglet.media.AudioFormat(num_channels, num_bits, rate)

        super(ArraySource, self).__init__(data.tostring(), audio_format)

    def _get_queue_source(self):
        return self

def playsound(data, rate):
    sound = ArraySource(rate, data)
    player = pyglet.media.Player()
    player.queue(sound)
    player.play()

def create_dtmf(f1, f2, fs, duration):
    t = np.linspace(0, duration, duration*fs)
    f = (np.sin(f1*np.pi*t)+np.sin(f2*np.pi*t))/2*scale
    return np.reshape(f, (len(f),1))

def encoder(symbol):
	for i in range(16):
		if symbol == keys[i]:
			f1 = F1[i/4] #row
			f2 = F2[i%4] #column
	data = range(FR)
	for i in range(FR):
		p = i*1.0/FR
		data[i]=int(scale+(np.sin(p*f1*2*np.pi)+np.sin(p*f2*2*np.pi))/2*scale)
	return data

if __name__ == '__main__':
    data = np.array(encoder('A'))
    data = np.reshape(data, (len(data),1))
    playsound(data, FR)
    plt.plot(data)
    plt.show()