#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : reflectivity_profile.py
Purpose : make reflectivity profile for erzsol3
Creation Date : 21-03-2017
Last Modified : Tue 21 Mar 2017 07:29:26 PM EDT
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
    write_reflectivity(reflect,'test')
    return reflect

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
    z = np.log(mantle[:,0]/r0)*r0
    rho = mantle[:,1]*r0/mantle[:,0]
    vp = mantle[:,2]*r0/mantle[:,0]
    vs = mantle[:,3]*r0/mantle[:,0]
    qp = mantle[:,4]
    qs = mantle[:,5]
    reflect = np.hstack((np.array([vp]).T,
              np.array([vs]).T,
              np.array([rho]).T,
              z.max()-np.array([z]).T,
              np.array([qp]).T,
              np.array([qs]).T,
              ))
    return reflect

def write_reflectivity(reflect,name):
    out = np.hstack((np.array([np.ones(reflect.shape[0])]).T,reflect))
    print reflect.shape

    with open(name+'.mod','w') as f:
        f.write(name+'\n')
        f.write('    {}       {}\n'.format(reflect.shape[0],str(1)))
        np.savetxt(f,out,delimiter="   ",fmt='%.1f')

reflect = main()
