#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : find_bounce_points.py
Purpose : find all reflection points for first order reverbs
Creation Date : 14-06-2017
Last Modified : Thu 15 Jun 2017 10:17:17 AM EDT
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
import geopy
from geopy.distance import VincentyDistance as vd
from scipy.signal import argrelextrema
from obspy.taup import TauPyModel
model = TauPyModel(model='prem')

def main():
    dirname = '/home/samhaug/work1/ScS_reverb_data/obspyDMT/mag_7.0_8.0/20030620_061938.a/processed/'
    f = h5py.File(dirname+'bounce_points.h5','w')
    st = obspy.read(dirname+'st_T.pk')
    st = seispy.filter.range_filter(st,(10,90))
    fam_list = make_families()
    dir_list = ['top','bot','top','bot','top']
    fam_bounce = []
    depth = 670
    for jdx,jj in enumerate(fam_list):
        for idx,ii in enumerate(jj):
            phase_list = [ii.replace('#',str(depth))]
            c_list = np.array(bounce_family(st,phase_list,dir_list[jdx],depth))
            f.create_dataset('fam_'+str(jdx)+'/'+phase_list[0],data=np.array(c_list))

def make_families():
    first_family = ['sSv#SScS',
                    'sScSSv#S']
    second_family = ['ScS^#ScS']
    third_family = ['sSv#SScSScS',
                    'sScSSv#SScS',
                    'sScSScSSv#S']
    fourth_family = ['ScS^#ScSScS',
                      'ScSScS^#ScS']
    fifth_family = ['sSv#SScSScSScS',
                    'sScSSv#SScSScS',
                    'sScSScSSv#SScS',
                    'sScSScSScSSv#S']
    return [first_family,second_family,third_family,fourth_family,fifth_family]

def bounce_family(st,phase_list,direction,depth):
    c_list = []
    for tr in st:
        c_list.append(_find_bounce_coord(tr,phase_list,direction,depth))
    return c_list

def _find_bounce_coord(tr,phase_list,direction,depth):
    arr = model.get_pierce_points_geo(source_depth_in_km=tr.stats.sac['evdp'],
                                      source_latitude_in_deg=tr.stats.sac['evla'],
                                      source_longitude_in_deg=tr.stats.sac['evlo'],
                                      receiver_latitude_in_deg=tr.stats.sac['stla'],
                                      receiver_longitude_in_deg=tr.stats.sac['stlo'],
                                      phase_list=phase_list)

    a = np.array([[jj for jj in ii] for ii in arr[0].pierce])
    if direction == 'top':
        ex = argrelextrema(a[:-3],np.greater)[0]
        for ii in ex:
            if a[ii,-3] != depth:
                continue
            else:
               return a[ii,-2],a[ii,-1]

    if direction == 'bot':
        ex = argrelextrema(a[:-3],np.less)[0]
        for ii in ex:
            if a[ii,-3] != depth:
                continue
            else:
               return a[ii,-2],a[ii,-1]

main()




