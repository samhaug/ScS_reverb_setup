#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : deep_zeroth_order.py
Purpose : Get zeroth order reverbs for deep phases by subtracting perturbed
          first order traces in a wavelet_glossary
Creation Date : 15-03-2017
Last Modified : Wed 15 Mar 2017 03:22:39 PM EDT
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

def main():
    wave_glos = '/home/samhaug/work1/ScS_reverb_sims/wave_glossary/'
    lkup = '/home/samhaug/work1/ScS_reverb_sims/lookup_tables/'
    st = obspy.read('/home/samhaug/work1/ScS_reverb_sims/mineos/prem_568_FJ/st_clean.pk')

    lkup_table = lkup+'NA_prem_568_20160130.h5'

    dict_670 = get_waveforms(wave_glos+'prem_568_FJ_20160130.h5')
    dict_410,dict_220 = scale_waveforms(dict_670,0.8,0.7)
    t_670,t_410,t_220 = find_times(lkup_table,dict_670,dict_410,dict_220)
    place_reverbs(st,dict_670,dict_410,dict_220,t_670,t_410,t_220)

def get_waveforms(h5):
    f = h5py.File(h5,'r')
    dict_670 = {}
    for keys in f:
        dict_670[keys] = f[keys][...]
    return dict_670

def scale_waveforms(dict_670,sc410,sc220):
    dict_410 = {}
    dict_220 = {}
    for keys in dict_670:
        dict_410[keys] = dict_670[keys]*sc410
        dict_220[keys] = dict_670[keys]*sc220
    return dict_410,dict_220

def find_times(lkup,dict_670,dict_410,dict_220,**kwargs):
    stat = kwargs.get('stat','DSXP/')

    f = h5py.File(lkup)
    ScS2 = f[stat+'ScS2'][...][0][1]

    t_670 = {}
    for keys in dict_670:
        i = np.argmin(np.abs(f[stat+keys][:,0]-670))
        t = f[stat+keys][i,1]
        t_670[keys] = int(t-ScS2)
    t_410 = {}
    for keys in dict_410:
        i = np.argmin(np.abs(f[stat+keys][:,0]-400))
        t = f[stat+keys][i,1]
        t_410[keys] = int(t-ScS2)
    t_220 = {}
    for keys in dict_220:
        i = np.argmin(np.abs(f[stat+keys][:,0]-220))
        t = f[stat+keys][i,1]
        t_220[keys] = int(t-ScS2)

    return t_670,t_410,t_220

def place_reverbs(st,dict_670,dict_410,dict_220,t_670,t_410,t_220):
    blank_410 = np.zeros(len(st[0].data))
    for keys in dict_410:
        t = 400+t_410[keys]
        r = dict_410[keys]
        blank_410[t:t+len(r)] = r

    blank_670 = np.zeros(len(st[0].data))
    for keys in dict_670:
        t = 400+t_670[keys]
        r = dict_670[keys]
        try:
            blank_670[t:t+len(r)] = r
        except ValueError:
            continue

    blank_220 = np.zeros(len(st[0].data))
    for keys in dict_220:
        t = 400+t_220[keys]
        r = dict_220[keys]
        blank_220[t:t+len(r)] = r

    plt.plot(blank_410+blank_220+blank_670,alpha=0.5)
    plt.plot(st[0].data,alpha=0.5)
    plt.show()

#def subtract_waveforms():
#def pass_zero_order():

main()
