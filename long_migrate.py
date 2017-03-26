#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : long_migrate.py
Purpose : Use respnse_array from h5file made with make_long_migrate_wave.py
          to create a depth migrated reflection profile for data
Creation Date : 20-03-2017
Last Modified : Fri 24 Mar 2017 03:03:45 PM EDT
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
    tr = st[10]
    response,depths = read_h5('test1.h5')
    R,depths = migrate(tr,response,depths)
    save_migrate(tr,R,depths)

def stream_setup():
    sim_dir = '/home/samhaug/work1/ScS_reverb_sims/mineos/'
    st = obspy.read(sim_dir+'prem_013016E/st_T.pk')
    st.integrate().detrend().integrate().detrend()
    st.interpolate(1)
    st.filter('bandpass',freqmin=1./85,freqmax=1./15,zerophase=True)
    st.normalize()
    for idx,tr in enumerate(st):
       st[idx] = seispy.data.phase_window(tr,phase=['ScSScS'],window=(-400,2400))
       st[idx].stats.sac['o'] += -1468
    st.sort(['location'])
    return st

def read_h5(name):
    f = h5py.File('/home/samhaug/work1/ScS_reverb_sims/long_migrate/'+name,'r')
    response = f['response'][...]
    depths = f['depths'][...]
    f.close()
    return response,depths

def migrate(tr,response,depths):
    data = tr.data
    R = []
    for ii in response:
        denom = np.dot(ii,ii)
        corr_sig = correlate(data,ii,mode='same')
        R.append(corr_sig[int(len(corr_sig)/2.)]/denom)
    plt.plot(R,depths)
    plt.ylim(depths.max(),depths.min())
    plt.axhline(670)
    plt.axhline(400)
    plt.axhline(220)
    plt.show()
    return R,depths

def save_migrate(tr,R,depths):
    print np.array([R]).T.shape,np.array([depths]).T.shape
    r = np.hstack((np.array([R]).T,np.array([depths]).T))
    np.savetxt('/home/samhaug/work1/ScS_reverb_sims/reflectivity_profiles/{}_{}.dat'.format(
               tr.stats.station,tr.stats.network),r)

main()
