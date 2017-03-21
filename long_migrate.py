#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : long_migrate.py
Purpose : Use respnse_array from h5file made with make_long_migrate_wave.py
          to create a depth migrated reflection profile for data
Creation Date : 20-03-2017
Last Modified : Tue 21 Mar 2017 10:49:45 AM EDT
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
from scipy.signal import correlate

def main():
    st = stream_setup()
    response,depths = read_h5('test1.h5')
    migrate(st[10].data,response,depths)

def stream_setup():
    sim_dir = '/home/samhaug/work1/ScS_reverb_sims/mineos/'
    st = obspy.read(sim_dir+'prem_468_FJ/st_T.pk')
    st.integrate().detrend().integrate().detrend()
    st.interpolate(1)
    st.filter('bandpass',freqmin=1./60,freqmax=1./10,zerophase=True)
    st.normalize()
    for idx,tr in enumerate(st):
       st[idx] = seispy.data.phase_window(tr,phase=['ScSScS'],window=(-400,2400))
       st[idx].stats.sac['o'] += -1468
    return st

def read_h5(name):
    f = h5py.File('/home/samhaug/work1/ScS_reverb_sims/long_migrate/'+name,'r')
    response = f['response'][...]
    depths = f['depths'][...]
    f.close()
    return response,depths

def migrate(data,response,depths):
    R = []
    for ii in response:
        corr_sig = correlate(data,ii,mode='same')
        R.append(corr_sig[int(len(corr_sig)/2.)])
    plt.plot(R,depths)
    plt.ylim(depths.max(),depths.min())
    plt.axhline(670)
    plt.axhline(400)
    plt.axhline(220)
    plt.show()

main()
