#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : long_migrate.py
Purpose : Use respnse_array from h5file made with make_long_migrate_wave.py
          to create a depth migrated reflection profile for data
Creation Date : 20-03-2017
Last Modified : Mon 03 Apr 2017 06:16:45 PM EDT
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
from scipy.signal import tukey
from obspy.taup import TauPyModel

def main():
    st = stream_setup('prem_12.0_10','VEL_PERT_RUNS/japan_111702/japan_111702_12.0/')
    tr = mask_zeroth_order('prem_12.0_10',st[1])
    response,depths = read_h5('test1.h5')
    R,depths = migrate(tr,response,depths)
    save_migrate(tr,R,depths)

def stream_setup(model,sim):
    model = TauPyModel(model=model)
    sim_dir = '/home/samhaug/work1/ScS_reverb_sims/mineos/'+sim
    st = obspy.read(sim_dir+'st_T.pk')
    st.integrate().detrend().integrate().detrend()
    st.interpolate(1)
    st.filter('bandpass',freqmin=1./75,freqmax=1./15,zerophase=True)
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

def mask_zeroth_order(model,tr):
    model = TauPyModel(model=model)
    phase_list = ['ScSScS',
                  'sScSScS',
                  'sScSScSScS',
                  'ScSScSScS',
                  'ScSScSScSScS',
                  'sScSScSScSScS']
    mlen = 50
    shift = -10
    arrival = model.get_travel_times(source_depth_in_km=tr.stats.sac['evdp'],
                                         distance_in_degree=tr.stats.sac['gcarc'],
                                         phase_list=phase_list)
    mask = np.ones(len(tr.data))
    for ii in arrival:
        t = int(400+ii.time-arrival[0].time)
        #print t
        try:
            mask[t+shift:t+mlen+shift] = 1-tukey(mlen,0.8)
        except ValueError:
            continue
    tr.data *= mask
    #plt.plot(tr.data)
    #plt.plot(tr.data*mask)
    #plt.show()
    return tr

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
