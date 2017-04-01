#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : long_migrate.py
Purpose : Use respnse_array from h5file made with make_long_migrate_wave.py
          to create a depth migrated reflection profile for data
Creation Date : 20-03-2017
Last Modified : Fri 31 Mar 2017 01:01:04 PM EDT
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
from obspy.taup import TauPyModel

def main():
    st = stream_setup('prem_9.0_10')
    tr = st[-7]
    response,depths = read_h5('test1.h5')
    R,depths = migrate(tr,response,depths)
    save_migrate(tr,R,depths)

def stream_setup(model):
    model = TauPyModel(model=model)
    sim_dir = '/home/samhaug/work1/ScS_reverb_data/'
    st = obspy.read(sim_dir+'013016_japan/st_T.pk')
    st.integrate().detrend()
    st.interpolate(1)
    st.filter('bandpass',freqmin=1./75,freqmax=1./10,zerophase=True)
    st.normalize()
    for idx,tr in enumerate(st):
       st[idx] = seispy.data.phase_window(tr,phase=['ScSScS'],window=(-400,2400))
       arrival = model.get_travel_times(source_depth_in_km=st[idx].stats.sac['evdp'],
                                        distance_in_degree=st[idx].stats.sac['gcarc'],
                                        phase_list=['ScSScS'])
       o = arrival[0].time-400.
       st[idx] = seispy.data.phase_window(tr,phase=['ScSScS'],window=(-400,2400))
       st[idx].stats.sac['o'] += -1*o
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
