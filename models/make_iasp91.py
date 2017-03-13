#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : make_iasp91.py
Purpose : Make iasp91 input file for mineos to use. Needs prem_noocean.txt
          as well as a csv IASP91 file from IRIS
Creation Date : 13-03-2017
Last Modified : Mon 13 Mar 2017 06:26:07 PM EDT
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
import scipy

def main():
    meta,cores,mantle = read_prem('prem_noocean.txt')
    irad,ivp,ivs = read_iasp91('IASP91.csv')
    imid,iend,irho = interp_iasp(irad,mantle[:,0],mantle[:,4:6],mantle[:,-1],mantle[:,1])
    imantle = np.hstack((
              np.array([irad]).T,
              np.array([irho]).T,
              np.array([ivp]).T,
              np.array([ivs]).T,
              imid,
              np.array([ivp]).T,
              np.array([ivs]).T,
              np.array([iend]).T
              ))

    inarray = np.vstack((cores,np.flipud(imantle)))
    write_model(inarray,'iasp91.txt')

def read_prem(infile):
    b = np.genfromtxt(infile,skip_header=3)
    rad = b[:,0]
    rho = b[:,1]
    mid = b[:,4:6]
    a = open(infile)
    f = a.read()
    meta = f.strip().split('\n')[2].split()
    knots,nic,noc = meta[0],meta[1],meta[2]
    a.close()
    model = np.genfromtxt(infile,skip_header=3)
    cores = model[0:int(noc),:]
    mantle = model[int(noc)::,:]
    return meta,cores,mantle

def read_iasp91(infile):
    a = np.genfromtxt(infile,delimiter=',')
    mantle = a[0:137,:]
    mantle[-1,1] = 3480.
    irad = mantle[:,1]*1000
    ivp = mantle[:,2]*1000
    ivs = mantle[:,3]*1000
    return irad,ivp,ivs

def interp_iasp(irad,prad,pmid,pend,prho):
    fmid = scipy.interpolate.interp1d(prad,pmid,axis=0)
    fend = scipy.interpolate.interp1d(prad,pend)
    frho = scipy.interpolate.interp1d(prad,prho)
    imid = fmid(irad)
    iend = fend(irad)
    irho = frho(irad)
    return imid,iend,irho

def write_model(inarray,title):
    with open(title, 'w') as f:
        f.write('IASP91 with PREM att and rho\n')
        f.write('1    1.0000   1\n')
        f.write('{}   {}   {}\n'.format(str(inarray.shape[0]),'33','66'))
        np.savetxt(f, inarray, delimiter="   ",fmt='%.1f')

main()



