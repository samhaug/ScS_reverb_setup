#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : make_iasp91.py
Purpose : Make iasp91 input file for mineos to use. Needs prem_noocean.txt
          as well as a csv IASP91 file from IRIS
Creation Date : 13-03-2017
Last Modified : Mon 13 Mar 2017 02:28:27 PM EDT
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
    prad,prho,pmid,pend = read_prem('prem_noocean.txt')
    irad,ivp,ivs = read_iasp91('IASP91.csv')
    pvp,pvs = interp_iasp(prad,irad,ivp,ivs)
    inarray = np.hstack((
              np.array([prad]).T,
              np.array([prho]).T,
              np.array([pvp]).T,
              np.array([pvs]).T,
              pmid,
              np.array([pvp]).T,
              np.array([pvs]).T,
              np.array([pend]).T
              ))
    write_model(inarray,'iasp91.txt')

def read_prem(infile):
    a = np.genfromtxt(infile,skip_header=3)
    rad = a[:,0]
    rho = a[:,1]
    mid = a[:,4:6]
    end = a[:,-1]
    return rad,rho,mid,end

def read_iasp91(infile):
    a = np.genfromtxt(infile,delimiter=',')
    rad = a[:,1]*1000
    vp = a[:,2]*1000
    vs = a[:,3]*1000
    return rad[::-1],vp,vs

def interp_iasp(prad,irad,vp,vs):
    fp = scipy.interpolate.interp1d(irad,vp)
    fs = scipy.interpolate.interp1d(irad,vs)
    prem_vs = fs(prad)
    prem_vp = fp(prad)
    return prem_vp, prem_vs

def write_model(inarray,title):
    with open(title, 'w') as f:
        f.write('IASP91 with PREM att and rho\n')
        f.write('1    1.0000   1\n')
        f.write('{}   {}   {}\n'.format(str(inarray.shape[0]),'33','66'))
        np.savetxt(f, inarray, delimiter="   ",fmt='%.1f')

main()



