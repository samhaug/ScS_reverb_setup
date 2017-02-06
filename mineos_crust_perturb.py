#!/home/samhaug/anaconda2/python

from sys import argv
import numpy as np
from scipy.signal import blackman
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

'''
Perturb the thickness and velocity of the crust. To be used for
waveform inversion of zeroth order reverberations to constrain moho depth
and reflection coefficient.

argv[1]: path to mineos model txt file. Must be in tabular format.
'''

def perturb_thickness_crust(filename,thickness,vs_pert):
    model = np.genfromtxt(filename,skip_header=3)
    no_crust = model[0:157,:]
    crust = model[157::,:]
    new_crust = np.arange(6371000-(thickness*1000),6372000,1000)
    crust_param = np.tile(crust[-1,1::],(len(new_crust),1))
    final_crust = np.hstack((np.transpose([new_crust]),crust_param))

    ind = np.argmin(np.abs(final_crust[0,0]-model[:,0]))
    print ind
    final_model = np.vstack((model[0:ind,:],final_crust))
    return final_model

def write_output(out_array,discont,amp,**kwargs):
    title= kwargs.get('title','test.txt')
    with open(title, 'w') as f:
        f.write('Radius: {} m, Amplitude: {} %vs\n'.format(str(discont),
                str(amp)))
        f.write('1    1.0000   1\n')
        f.write('{}   {}   {}\n'.format(str(out_array.shape[0]),33,66))
        np.savetxt(f, out_array, delimiter="   ",fmt='%.1f')


a = perturb_thickness_crust('./models/prem_noocean.txt',10,None)
b = perturb_thickness_crust('./models/prem_noocean.txt',40,None)



