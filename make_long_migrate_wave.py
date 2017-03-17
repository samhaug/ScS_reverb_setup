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

def main():
    st = stream_setup()
    lkup = lkup_setup()
    #seispy.plot.plot(st[10],phase_list=['ScSScS','ScS^670ScSScS'])
    mask_660(st[10],lkup)

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

def mask_660(tr,lkup):
    evdp = tr.stats.sac['evdp']
    gcarc = tr.stats.sac['gcarc']
    model = TauPyModel('prem50')
    stat = tr.stats.station

    # t is for table
    t = lkup[stat]
    ScS2 = t['ScS2'][0,1]
    sScS2 = t['sScS2'][0,1]
    dp = sScS2-ScS2

    master_mask = []
    for keys in t:
        if 'b' in t or 't' in t:
            mask = np.zeros(len(tr.data))
            i = np.argmin(np.abs(t[keys][:,0]-670))
            # r to sr is the prem predicted depth phase separation for reverb
            r = t[keys][i,1]-ScS2+400
            sr = dp+r
            mask[r:sr] = 1.0
            master_mask.append(mask)
    print master_mask
    plt.plot(tr.data)
    for ii in master_mask:
        plt.plot(ii,alpha=0.5)
    plt.show()
    return np.sum(master_mask,axis=0)

#def shift_wavelets(tr,depth,lkup):
#def write_wave_glossary(h5file,key,data):

main()


