#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : check_travel_times.py
Purpose : make traveltime plot of first order ScS reverberations with respect
         to ScS_2. Make this as a function of earthquake depth to find the best
         depth range.
Creation Date : 26-03-2017
Last Modified :
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
model = TauPyModel(model='prem')

def main():
    d_list = np.linspace(100,670,num=50)
    phase_list = ['ScSScSScS',
                  'sScSScSScS',
                  'ScSScS^670ScS',
                  'sScSScS^670ScS',
                  'ScSScS^400ScS',
                  'sScSScS^400ScS',
                  'ScSScS^220ScS',
                  'sScSScS^220ScS']

    for ii in phase_list:
        z,n = make_lookup_list(d_list,ii)
        plt.plot(d_list,z,label=n)
    plt.legend()
    plt.show()


def lookup_times(d,phase):
    zero_arrivals = model.get_travel_times(source_depth_in_km=d,
                                      distance_in_degree=30,
                                      phase_list = [phase])
    return zero_arrivals

def make_lookup_list(d_list,phase):
    zero_time_list = []
    for d in d_list:
        zero_arrivals = lookup_times(d,phase)
        zero_depth = []
        for ii in zero_arrivals:
            zero_depth.append(ii.time)
        zero_time_list.append(zero_depth)

    return np.array(zero_time_list),zero_arrivals[0].name

d_list,zt,ft,zn,fn = main()




