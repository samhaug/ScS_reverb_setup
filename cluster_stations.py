#!/home/samhaug/anaconda2/bin/python 
'''
==============================================================================

File Name : cluster_stations.py
Purpose : cluster stations into regions within a circle
Creation Date : 31-05-2017
Last Modified : Wed 31 May 2017 12:00:02 PM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
from matplotlib import pyplot as plt
from subprocess import call
from os import listdir
import h5py
import obspy
from obspy.geodetics import gps2dist_azimuth as gps
import seispy
from sklearn.cluster import k_means

def main():
    #st = read_stream('/home/samhaug/work1/ScS_reverb_data/2015-11-24-mww76-peru-brazil-border-region-5/')
    st = read_stream('/home/samhaug/work1/ScS_reverb_data/031011G_sac/')
    st_out = cluster(st)
    write_stream(st_out,'/home/samhaug/work1/ScS_reverb_data/031011G_sac/')

def read_stream(dirname):
    print 'reading'
    st = obspy.read(dirname+'*t.sac')
    st.filter('bandpass',freqmin=1./80,freqmax=1/15)
    for idx,tr in enumerate(st):
        st[idx] = seispy.data.phase_window(tr,['ScSScS'],
                                           window=(-100,800))
    return st

def write_stream(st,dirname):
    print ('writing')
    st.write(dirname+'k_means.pk',format='PICKLE')

def cluster(st,arc=4,k_number=100):
    print 'clustering'
    st_out = st[0:k_number].copy()
    for tr in st_out:
        tr.stats.sac['data_list'] = []
    coord_arr = np.array([[tr.stats.sac['stla'],
                            tr.stats.sac['stlo']] for tr in st])
    K = k_means(coord_arr,k_number)

    cluster_list = K[0]
    for idx,ii in enumerate(cluster_list):
        st_out[idx].stats.station = 'C'+str(idx)
        st_out[idx].stats.network = 'KM'
        st_out[idx].stats.sac['stla'] = ii[0]
        st_out[idx].stats.sac['stlo'] = ii[1]
        for tr in st:
            a = gps(ii[0],ii[1],tr.stats.sac['stla'],
                                tr.stats.sac['stlo'])[0]/111195.
            if a < arc:
                st_out[idx].stats.sac['data_list'].append(tr.data)
            else:
                continue
        st_out[idx].data = np.mean(st_out[idx].stats.sac['data_list'],axis=0)
        st_out[idx].std = np.std(st_out[idx].stats.sac['data_list'],axis=0)
    return st_out

main()


