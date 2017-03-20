#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : make_full_migrate_wave.py
Purpose : Make full first order reverberation interval to be used in
          long_migrate_reverb. First mask all but the 660 reverb. Then add or
          take away whitespace inbetween wavelets to simulate different depths.
Creation Date : 17-03-2017
Last Modified : Fri 17 Mar 2017 11:37:01 AM EDT
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
from obspy.taup import TauPyModel
from scipy.signal import tukey
from scipy.signal import correlate

def main():
    st = stream_setup()
    lkup = lkup_setup()
    tr = st[10]
    master_mask = mask_670(tr,lkup,plot=False)
    new_mask = shift_depth(tr,master_mask,368)
    response_list = shift_discont(tr,new_mask,lkup)
    for idx,ii in enumerate(response_list[::10]):
        plt.plot(idx+ii/ii.max(),alpha=0.5,color='k')
    plt.show()

def lkup_setup():
    lkup_dir = '/home/samhaug/work1/ScS_reverb_sims/lookup_tables/'
    lkup = h5py.File(lkup_dir+'NA_prem_368_20160130.h5','r')
    return lkup

def stream_setup():
    sim_dir = '/home/samhaug/work1/ScS_reverb_sims/mineos/'
    st = obspy.read(sim_dir+'prem_368_FJ/st_T.pk')
    st.integrate().detrend().integrate().detrend()
    st.interpolate(1)
    st.filter('bandpass',freqmin=1./75,freqmax=1./15,zerophase=True)
    st.normalize()
    for idx,tr in enumerate(st):
       st[idx] = seispy.data.phase_window(tr,phase=['ScSScS'],window=(-400,2400))
       st[idx].stats.sac['o'] += -1468
    return st

def mask_670(tr,lkup,**kwargs):
    plot = kwargs.get('plot',False)
    evdp = tr.stats.sac['evdp']
    gcarc = tr.stats.sac['gcarc']
    stat = tr.stats.station
    #mlen is length of wavelet window.
    mlen = 35
    #mshi is offset time of window sampler.
    mshi = 5

    # t is for table
    t = lkup[stat]
    ScS2 = t['ScS2'][0,1]
    sScS2 = t['sScS2'][0,1]
    dp = sScS2-ScS2

    master_mask = {}
    if plot:
        plt.plot(tr.data,alpha=0.3,color='k')
    for keys in t:
        if keys[0] == 'b' or keys[0] == 't':
            mask = np.zeros(len(tr.data))
            dmask = np.zeros(len(tr.data))
            i = np.argmin(np.abs(t[keys][:,0]-670))
            #r to sr is the prem predicted depth phase separation for reverb
            r = t[keys][i,1]-ScS2+400
            sr = dp+r
            try:
                mask[r+mshi:r+mshi+mlen] += 1.0*tukey(mlen,0.20)
                dmask[sr+mshi:sr+mshi+mlen] += 1.0*tukey(mlen,0.20)
                master_mask[keys] = (mask*tr.data)
                master_mask['d'+keys] = (dmask*tr.data)
                if plot:
                    plt.plot(tr.data*dmask,color='b')
                    plt.plot(tr.data*mask,color='b')
            except ValueError:
                continue
    if plot:
        plt.show()
    return master_mask

def shift_depth(tr,master_mask,depth):
    model = TauPyModel('prem50')
    gcarc = tr.stats.sac['gcarc']
    evdp  = tr.stats.sac['evdp']
    a = model.get_travel_times(source_depth_in_km=evdp,distance_in_degree=gcarc,
                           phase_list=['ScSScS','sScSScS'])
    b = model.get_travel_times(source_depth_in_km=depth,distance_in_degree=gcarc,
                           phase_list=['ScSScS','sScSScS'])
    shift = int(b[1].time-b[0].time)-int(a[1].time-a[0].time)
    print 'shift time:', shift
    new_mask = {}
    for keys in master_mask:
        if keys[0] == 'd':
            master_mask[keys]= seispy.data.roll_zero(master_mask[keys],shift)
    for keys in master_mask:
        if keys[0] == 'd':
            continue
        new_mask[keys] = master_mask[keys]+master_mask['d'+keys]
    return new_mask

def shift_discont(tr,new_mask,lkup):
    stat = tr.stats.station
    depths = np.arange(10,2800,2)
    response_list = []

    for d in depths:
        single_response = []
        for lkeys in lkup[stat]:
            if lkeys == 'ScS2' or lkeys == 'sScS2':
                continue
            i = np.argmin(np.abs(lkup[stat+'/'+lkeys][:,0]-d))
            td = lkup[stat+'/'+lkeys][i,1]
            i670 = np.argmin(np.abs(lkup[stat+'/'+lkeys][:,0]-670))
            t670 = lkup[stat+'/'+lkeys][i670,1]
            shift = int(td-t670)
            try:
                single_response.append(seispy.data.roll_zero(new_mask[lkeys],shift))
            except KeyError:
                continue
        response_list.append(np.sum(single_response,axis=0))
    response_array = np.array(response_list)
    return response_array

main()
