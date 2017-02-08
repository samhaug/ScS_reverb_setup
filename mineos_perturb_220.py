#!/home/samhaug/anaconda2/python

from sys import argv
import numpy as np
from scipy.signal import blackman
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

'''
Use this to perturb the upper 220km of the mantle to get best estimate
of ScS reverberation times.

This script wants to read prem_noocean.txt
'''

def read_model(filename,perturb):
    p = np.genfromtxt(filename,skip_header=3)
    p[147::,3] *= (1+perturb/100.)
    p[147::,7] *= (1+perturb/100.)
    return p

def write_output(out_array,amp):
    #title = kwargs.get('title','test.txt')
    title = './perturb_220/'+'220_'+str(amp)+'.txt'
    with open(title, 'w') as f:
        f.write('Radius: {} km, Amplitude: {} %vs\n'.format(str(220),
                str(amp)))
        f.write('1    1.0000   1\n')
        f.write('{}   {}   {}\n'.format(str(out_array.shape[0]),33,66))
        np.savetxt(f, out_array, delimiter="   ",fmt='%.1f')

amp = 0.4
out_array = read_model('./models/prem_noocean.txt',amp)
write_output(out_array,amp)




