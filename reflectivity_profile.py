#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : reflectivity_profile.py
Purpose : make reflectivity profile for erzsol3
Creation Date : 21-03-2017
Last Modified : Wed 22 Mar 2017 11:11:02 AM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
from matplotlib import pyplot as plt
from subprocess import call
from os import listdir
import h5py
import obspy
import seispy

def main():
    meta,cores,mantle = read_model('./models/noocean_notrans.txt')
    reflect = flat_earth_approx(mantle)
    thick = depth_to_thick(reflect)
    write_reflectivity(thick,'test')
    #return out

def read_model(filename):
    a = open(filename)
    f = a.read()
    meta = f.strip().split('\n')[2].split()
    knots,nic,noc = meta[0],meta[1],meta[2]
    print noc
    a.close()
    model = np.genfromtxt(filename,skip_header=3)
    cores = model[0:int(noc),:]
    mantle = model[int(noc)::,:]
    return meta,cores,mantle

def flat_earth_approx(mantle):
    r0 = mantle[:,0].max()
    z = (np.log(mantle[:,0]/r0)*r0)/1000.
    rho = (mantle[:,1]*r0/mantle[:,0])/1000.
    vp = (mantle[:,2]*r0/mantle[:,0])/1000.
    vs = (mantle[:,3]*r0/mantle[:,0])/1000.
    qp = 1./mantle[:,4]
    qs = 1./mantle[:,5]
    reflect = np.hstack((np.array([vp]).T,
              np.array([vs]).T,
              np.array([rho]).T,
              z.max()-np.array([z]).T,
              np.array([qp]).T,
              np.array([qs]).T,
              ))
    return reflect

def depth_to_thick(reflect):
    reflect = np.flipud(reflect)
    thick = reflect.copy()
    for idx,ii in enumerate(reflect):
        try:
            thick[idx][3] = reflect[idx+1][3]-reflect[idx][3]
        except IndexError:
            continue
    return thick[0:-1]

def write_reflectivity(reflect,name):
    out = np.hstack((np.array([np.zeros(reflect.shape[0])]).T,reflect))

    with open(name+'.mod','w') as f:
        f.write(name+'\n')
        f.write('    {}       {}\n'.format(reflect.shape[0],str(1)))
        #np.savetxt(f,out,delimiter="   ",fmt='%.5f')
        for ii in out:
            f.write('{0:.0f}   {1:.3f}   {2:.3f}   {3:.3f}   {4:.3f}   {5:.3f}   {6:.3f}\n'.format(
                    ii[0],ii[1],ii[2],ii[3],ii[4],ii[5],ii[6]))
        #np.savetxt(f,out,delimiter="   ")


main()






