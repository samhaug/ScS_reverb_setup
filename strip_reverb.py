#!/home/samhaug/anaconda2/bin/python

import obspy
import numpy as np
from scipy.signal import tukey
from obspy.taup import TauPyModel
from matplotlib import pyplot as plt

model = TauPyModel(model='prem')

'''
This script will strip the 660 reverberations from an AxiSEM synthetic.
'''

def find_ttimes(tr):
    '''
    Find the traveltimes of ScS reverberations for the 660
    '''
    print('find_ttimes')
    gcarc = tr.stats.sac['gcarc']
    evdp = tr.stats.sac['evdp']
    time_d = {}
    phase_list = [
                  'ScS^660ScS',
                  'sScS^660ScS',

                  'ScSSv660SScS',
                  'sScSSv660SScS',

                  'ScS^660ScSScS',
                  'sScS^660ScSScS',

                  'ScSSv660SScSScS',
                  'sScSSv660SScSScS',

                  'ScS^660ScSScSScS',
                  'sScS^660ScSScSScS',

                  'ScSSv660SScSScSScS',
                  'sScSSv660SScSScSScS',

                  'ScS^660ScSScSScSScS',
                  'sScS^660ScSScSScSScS'
                 ]

    arrivals = model.get_travel_times(source_depth_in_km=evdp,
                                   distance_in_degree=gcarc,
                                   phase_list = phase_list)
    for i,j in enumerate(phase_list):
        time_d[j] = arrivals[i].time

    return time_d

def strip_reverb(tr,time_d):
    '''
    Use a tukey window to mask an axisem
    trace and only pass the reverberations
    '''
    sr = tr.stats.sampling_rate
    npts = tr.stats.npts

    def make_mask(phase):
        S = time_d[phase]
        sS = time_d['s'+phase]
        diff = int(abs(sS-S))
        window = (S-30,S+2*diff+30)
        t = tukey(int((window[1]-window[0])*sr),0.2)
        start = np.zeros(int(window[0]*sr))
        end = np.zeros(npts-int(window[1]*sr))
        mask = np.hstack((start,t,end))
        #plt.plot(tr.data*mask)
        #plt.show()
        if len(mask) > len(tr.data):
            mask = mask[0:len(mask)-len(tr.data)]
        if len(mask) < len(tr.data):
            mask = np.hstack((mask,np.zeros(len(tr.data)-len(mask))))
        return mask

    m1 = make_mask('ScS^660ScS')
    m2 = make_mask('ScS^660ScSScS')
    m3 = make_mask('ScS^660ScSScSScS')
    m4 = make_mask('ScSSv660SScSScS')
    m5 = make_mask('ScSSv660SScS')
    m6 = make_mask('ScS^660ScSScSScSScS')
    m7 = make_mask('ScSSv660SScSScSScS')
    plt.plot((m1+m2+m3+m4+m5+m6+m7)*tr.data,alpha=0.5)
    plt.plot(tr.data,alpha=0.5)
    plt.show()


st = obspy.read('/home/samhaug/work1/ScS_reverb_sims/axisem/FJ_prem/st_T.pk')
st.filter('lowpass',freq=1./25)
st.normalize()

time_d = find_ttimes(st[10])
strip_reverb(st[10],time_d)



