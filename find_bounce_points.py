#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : find_bounce_points.py
Purpose : find all reflection points for first order reverbs
Creation Date : 14-06-2017
Last Modified : Wed 14 Jun 2017 06:24:29 PM EDT
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
    st = obspy.read('/home/samhaug/work1/ScS_reverb_data/obspyDMT/mag_7.0_8.0/20030620_061938.a/processed/st_T.pk')
    fam1,fam2,fam3,fam4,fam5 = make_families()
    fam_bounce = []
    color = ['b','r','g','k']
    for idx,ii in enumerate(fam3):
        phase_list = [ii.replace('#','670')]
        c_list = np.array(bounce_family(st,phase_list,670))
        plt.scatter(c_list[:,1],c_list[:,0],marker='+',color=color[idx])
    plt.show()

def make_families():
    first_family = ['sSv#SScS',
                    'sScSSv#S']
    second_family = ['ScS^#ScS']
    third_family = ['sSv#SScSScS',
                    'sScSv#SScS',
                    'sScSScSSv#S']
    fourth_family = ['ScSScSScS',
                     'ScS^#ScSScS',
                     'ScSScS^#ScS']
    fifth_family = ['sSv#SScSScSScS',
                    'sScSv#SScSScS',
                    'sScSScSv#SScS',
                    'sScSScSScSSv#S']
    return first_family,second_family,third_family,fourth_family,fifth_family

def bounce_family(st,phase_list,depth):
    c_list = []
    for tr in st:
        c_list.append(_find_bounce_coord(tr,phase_list,'top',depth))
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




