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
from scipy.signal import cosine
from scipy.signal import tukey
from scipy.signal import correlate
from scipy.signal import fftconvolve
from scipy.signal import gaussian
from scipy.optimize import curve_fit

def main():
    model = 'prem_9.0_10'
    shift = 585
    st = stream_setup(model,'VEL_PERT_RUNS/japan_100113/japan_100113_9.0/st_T.pk')
    lkup = lkup_setup('japan_100113_v9.0.h5')
    tr = st[1]

    #master_mask,gauss_fit_mask = mask_670(tr,lkup,plot=True)
    master_mask,gauss_fit_mask = mask_670_single(tr,lkup,plot=True)
    switch_mask = shift_depth(tr,master_mask,shift,model)
    #comment out this line unless switching bottom with topside
    #switch_mask = t_b_switch(switch_mask)

    response_array,depths = shift_discont(tr,switch_mask,lkup)
    write_h5(response_array,depths,'test1.h5')

    for idx,ii in enumerate(response_array[::10]):
        plt.plot(idx+ii/ii.max(),alpha=0.5,color='k')
    plt.show()

def write_h5(response_array,depths,name):
    f = h5py.File('/home/samhaug/work1/ScS_reverb_sims/long_migrate/'+name,'w')
    f.create_dataset('response',data=response_array)
    f.create_dataset('depths',data=depths)
    f.close()

def lkup_setup(lkup_table):
    lkup_dir = '/home/samhaug/work1/ScS_reverb_sims/lookup_tables/'
    lkup = h5py.File(lkup_dir+lkup_table,'r')
    return lkup

def stream_setup(model,stream):
    model = TauPyModel(model=model)
    sim_dir = '/home/samhaug/work1/ScS_reverb_sims/mineos/'
    st = obspy.read(sim_dir+stream)
    st.integrate().detrend().integrate().detrend()
    st.interpolate(1)
    st.filter('bandpass',freqmin=1./80,freqmax=1./10,zerophase=True)
    st.normalize()
    for idx,tr in enumerate(st):
       arrival = model.get_travel_times(source_depth_in_km=st[idx].stats.sac['evdp'],
                                        distance_in_degree=st[idx].stats.sac['gcarc'],
                                        phase_list=['ScSScS'])
       o = arrival[0].time-400.
       st[idx] = seispy.data.phase_window(tr,phase=['ScSScS'],window=(-400,2400))
       st[idx].stats.sac['o'] += -1*o
    #seispy.plot.plot(st[10],phase_list=['ScSScS','ScS^670ScSScS'])
    #seispy.plot.plot(st[20],phase_list=['ScSScS','ScS^670ScSScS'])
    #seispy.plot.plot(st[0],phase_list=['ScSScS','ScS^670ScSScS'])
    st.sort(['location'])
    return st

def mask_670_single(tr,lkup,**kwargs):
    plot = kwargs.get('plot',False)
    evdp = tr.stats.sac['evdp']
    gcarc = tr.stats.sac['gcarc']
    stat = tr.stats.station
    #mlen is length of wavelet window.
    mlen = 90
    #mshi is offset time of window sampler.
    mshi = -30

    # t is for table
    t = lkup[stat]
    ScS2 = t['ScS2'][0,1]
    sScS2 = t['sScS2'][0,1]
    dp = sScS2-ScS2

    master_mask = {}
    gauss_fit_mask = {}
    if plot:
        plt.plot(tr.data,alpha=0.3,color='k')
    for keys in t:
        if keys[-1] == '4' and keys[0] == 't':
            continue
        if keys[0] == 'b':
            mask = np.zeros(len(tr.data))
            dmask = np.zeros(len(tr.data))
            i = np.argmin(np.abs(t[keys][:,0]-670))
            #r to sr is the prem predicted depth phase separation for reverb
            r = t[keys][i,1]-ScS2+400
            sr = dp+r
            try:
                #mask[int(r+mshi):int(r+mshi+mlen)] += 1.0*cosine(int(mlen))**2
                mask[int(r+mshi):int(r+mshi+mlen)] += 1.0*tukey(int(mlen),0.6)
                dmask[int(sr+mshi):int(sr+mshi+mlen)] += np.zeros(int(mlen))
                master_mask[keys] = (mask*tr.data)
                master_mask['d'+keys] = (dmask*tr.data)
                gauss_fit_mask[keys] = (mask*tr.data)[int(r+mshi):int(r+mshi+mlen)]
                gauss_fit_mask['d'+keys] = (mask*tr.data)[int(sr+mshi):int(sr+mshi+mlen)]
                if plot:
                    plt.plot(dmask*0.01,color='r')
                    plt.plot(mask*0.01,color='r')
                    plt.plot(tr.data*dmask,color='b')
                    plt.plot(tr.data*mask,color='b')
            except ValueError:
                continue
        if keys[0] == 't':
            mask = np.zeros(len(tr.data))
            dmask = np.zeros(len(tr.data))
            i = np.argmin(np.abs(t[keys][:,0]-670))
            #r to sr is the prem predicted depth phase separation for reverb
            r = int(t[keys][i,1]-ScS2+400)
            sr = int(dp+r)
            try:
                mask[r+mshi:r+mshi+mlen] += np.zeros(int(mlen))
                #dmask[sr+mshi:sr+mshi+mlen] += 1.0*cosine(int(mlen))**2
                dmask[sr+mshi:sr+mshi+mlen] += 1.0*tukey(int(mlen),0.6)
                master_mask[keys] = (mask*tr.data)
                master_mask['d'+keys] = (dmask*tr.data)
                gauss_fit_mask[keys] = (mask*tr.data)[r+mshi:r+mshi+mlen]
                gauss_fit_mask['d'+keys] = (mask*tr.data)[sr+mshi:sr+mshi+mlen]
                if plot:
                    plt.plot(dmask*0.01,color='r')
                    plt.plot(mask*0.01,color='r')
                    plt.plot(tr.data*dmask,color='b')
                    plt.plot(tr.data*mask,color='b')
            except ValueError:
                continue

    if plot:
        plt.tight_layout()
        plt.show()
    return master_mask,gauss_fit_mask

def mask_670(tr,lkup,**kwargs):
    plot = kwargs.get('plot',False)
    evdp = tr.stats.sac['evdp']
    gcarc = tr.stats.sac['gcarc']
    stat = tr.stats.station
    #mlen is length of wavelet window.
    mlen = 30
    #mshi is offset time of window sampler.
    mshi = 0

    # t is for table
    t = lkup[stat]
    ScS2 = t['ScS2'][0,1]
    sScS2 = t['sScS2'][0,1]
    dp = sScS2-ScS2

    master_mask = {}
    gauss_fit_mask = {}
    if plot:
        plt.plot(tr.data,alpha=0.3,color='k')
    for keys in t:
        if keys[-1] == '4' and keys[0] == 't':
            continue
        if keys[0] == 'b' or keys[0] == 't':
            mask = np.zeros(len(tr.data))
            dmask = np.zeros(len(tr.data))
            i = np.argmin(np.abs(t[keys][:,0]-670))
            #r to sr is the prem predicted depth phase separation for reverb
            r = t[keys][i,1]-ScS2+400
            sr = dp+r
            try:
                #mask[r+mshi:r+mshi+mlen] += 1.0*cosine(mlen)**2
                #dmask[sr+mshi:sr+mshi+mlen] += 1.0*cosine(mlen)**2
                mask[r+mshi:r+mshi+mlen] += 1.0*tukey(mlen,0.4)
                dmask[sr+mshi:sr+mshi+mlen] += 1.0*tukey(mlen,0.4)
                master_mask[keys] = (mask*tr.data)
                master_mask['d'+keys] = (dmask*tr.data)
                gauss_fit_mask[keys] = (mask*tr.data)[r+mshi:r+mshi+mlen]
                gauss_fit_mask['d'+keys] = (mask*tr.data)[sr+mshi:sr+mshi+mlen]
                if plot:
                    plt.plot(dmask*0.01,color='r')
                    plt.plot(mask*0.01,color='r')
                    plt.plot(tr.data*dmask,color='b')
                    plt.plot(tr.data*mask,color='b')
            except ValueError:
                continue
    if plot:
        plt.tight_layout()
        plt.show()
    return master_mask,gauss_fit_mask

def t_b_switch(new_mask):
    switch_mask = new_mask.copy()
    for keys in new_mask:
        if keys[0] == 't':
            b_key = 'b'+keys[1::]
            b = -1*new_mask[b_key]
            t = new_mask[keys]
            ratio = np.abs(t).max()/np.abs(b).max()
            cor = fftconvolve(t,b[::-1])
            midpoint = cor.shape[0]/2
            imax = np.where(cor == cor.max())[0][0]
            roll = -1*(midpoint-imax)
            t = np.roll(b,roll)
            switch_mask[keys] = t
        if keys[0] == 'b':
            switch_mask[keys] = new_mask[keys]
            continue
    return switch_mask

def shift_depth(tr,master_mask,depth,model):
    model = TauPyModel(model=model)
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
    depths = np.arange(30,2800,2)
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
    return response_array,depths

main()



