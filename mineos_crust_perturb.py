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

Must reference prem_flattop.txt

argv[1]: path to mineos model txt file. Must be in tabular format.
'''

def perturb_thickness_crust(filename,thickness,vs_pert):
    model = np.genfromtxt(filename,skip_header=3)
    no_crust = model[0:154,:]
    crust = model[154::,:]
    f = interp1d(crust[:,0],crust,axis=0)
    new_rad = np.linspace(crust[:,0][0],crust[:,0][-1],40)
    new_crust = f(new_rad)
    crust_perturb = np.ones(new_crust.shape[0])

    ind = np.argmin(np.abs(6371000.-(thickness*1000)-new_crust[:,0]))
    repeat_layer = new_crust[ind].copy()
    repeat_layer[3] = new_crust[0,3]
    print repeat_layer
    crust_perturb[ind::] += vs_pert/100.
    new_crust[:,3]*=crust_perturb

    new_crust = np.vstack((new_crust[:ind,:],
                           repeat_layer,
                           new_crust[ind::,:]))

    new_crust = np.vstack((no_crust,new_crust))
    np.savetxt('new_crust.txt', new_crust, delimiter="   ",fmt='%.0f')
    return new_crust

def write_output(out_array,discont,amp,**kwargs):
    title= kwargs.get('title','test.txt')
    with open(title, 'w') as f:
        f.write('Radius: {} m, Amplitude: {} %vs\n'.format(str(discont),
                str(amp)))
        f.write('1    1.0000   1\n')
        f.write('{}   {}   {}\n'.format(str(out_array.shape[0]),33,66))
        np.savetxt(f, out_array, delimiter="   ",fmt='%.1f')

out = perturb_thickness_crust('./models/prem_flattop.txt',35,3)
write_output(out,35000,3,title='35_3.txt')






