#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : set_bounce_points.py
Purpose : find all reflection points for first order reverbs and set in tr.dict
Creation Date : 14-06-2017
Last Modified : Sat 24 Jun 2017 01:06:02 PM EDT
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
    dirname = '/home/samhaug/work1/ScS_reverb_data/2017-06-14-mww69-guatemala/'
    f = h5py.File(dirname+'bounce_points.h5','w')
    st = obspy.read(dirname+'st_T_clean.pk')
    st = seispy.filter.range_filter(st,(10,90))
    fam_list = make_families()
    dir_list = ['top','bot','top','bot','top']
    fam_bounce = []
    depth = 670
    for jdx,jj in enumerate(fam_list):
        for idx,ii in enumerate(jj):
            phase_list = [ii.replace('#',str(depth))]
            for tdx,tr in enumerate(st):
                st[idx] = find_bounce_coord(tr,'fam_'+str(jdx),phase_list,dir_list[jdx],depth)
    st.write(dirname+'st_T_bounce.pk',format='PICKLE')

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

def find_bounce_coord(tr,family,phase_list,direction,depth):
    arr = model.get_pierce_points_geo(source_depth_in_km=tr.stats.sac['evdp'],
                                      source_latitude_in_deg=tr.stats.sac['evla'],
                                      source_longitude_in_deg=tr.stats.sac['evlo'],
                                      receiver_latitude_in_deg=tr.stats.sac['stla'],
                                      receiver_longitude_in_deg=tr.stats.sac['stlo'],
                                      phase_list=phase_list)

    tr.stats.bd = {}
    tr.stats.bd['fam_0'] = {}
    tr.stats.bd['fam_1'] = {}
    tr.stats.bd['fam_2'] = {}
    tr.stats.bd['fam_3'] = {}
    tr.stats.bd['fam_4'] = {}

    a = np.array([[jj for jj in ii] for ii in arr[0].pierce])
    if direction == 'top':
        ex = argrelextrema(a[:-3],np.greater)[0]
        for ii in ex:
            if a[ii,-3] != depth:
                continue
            else:
               tr.stats.bd[family][phase_list[0]] = [a[ii,-2],a[ii,-1]]
               #return a[ii,-2],a[ii,-1]

    if direction == 'bot':
        ex = argrelextrema(a[:-3],np.less)[0]
        for ii in ex:
            if a[ii,-3] != depth:
                continue
            else:
               tr.stats.bd[family][phase_list[0]] = [a[ii,-2],a[ii,-1]]
               #return a[ii,-2],a[ii,-1]
    return tr

main()




