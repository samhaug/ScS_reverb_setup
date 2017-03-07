#!/home/samhaug/anaconda2/bin/python

import numpy as np
from matplotlib import pyplot as plt
import obspy
import seispy

st = obspy.read('/home/samhaug/work1/ScS_reverb_sims/mineos/FJ/FJ_0.0/st_T.pk')
st.integrate().detrend().integrate().detrend()
st.interpolate(1)
std = obspy.read('/home/samhaug/work1/ScS_reverb_data/20160130/FJ.pk')
std.integrate().detrend()
std.interpolate(1)

std.filter('bandpass',freqmin=1./70,freqmax=1./15,zerophase=True)
std = seispy.data.align_on_phase(std,phase=['ScSScS'])
for idx,tr in enumerate(std):
    std[idx] = seispy.data.phase_window(tr,phase=['ScSScS'],window=(-400,2400))
    std[idx].stats.starttime=0
    std[idx].data += -1*std[idx].data.mean()
std.normalize()

st.filter('bandpass',freqmin=1./70,freqmax=1./15,zerophase=True)
st = seispy.data.align_on_phase(st,phase=['ScSScS'])
for idx,tr in enumerate(st):
    st[idx] = seispy.data.phase_window(tr,phase=['ScSScS'],window=(-400,2400))
    st[idx].stats.starttime=0
    st[idx].data += -1*st[idx].data.mean()
st.normalize()


plt.show()

