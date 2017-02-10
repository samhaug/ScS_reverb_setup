#!/home/samhaug/anaconda2/bin/python

import numpy as np

'''
Once you find an optimal perturbed model for explaining ScS_n traveltimes,
use this to make a Taup tvel file to predict traveltimes
'''

def read_mineos_vel(file_name):

    vel = np.genfromtxt(file_name,skip_header=3)
    vel = vel[:,0:4]
    vel *= 1/1000.
    vel[:,0] = 6371. - vel[:,0]
    vel = np.flipud(vel)
    vel = np.vstack((vel[:,0],vel[:,2],vel[:,3],vel[:,1])).T
    print vel.shape
    return vel

def write_taup(vel_array,fname):

    with open(fname,'w') as f:
        f.write('prem_model.tvel P\n')
        f.write('prem_model.tvel S\n')
        np.savetxt(f,vel_array,delimiter='   ',fmt='%0.4f')

vel = read_mineos_vel('perturb_220/220_2.5.txt')
write_taup(vel,'prem_p2.5.tvel')
