#!/home/samhaug/anaconda2/bin/python

import obspy
from obspy.taup import TauPyModel
import numpy as np
from glob import glob

'''
Use this to make lookup table for ScS reverb traveltime perturbations.

'''

def reverb_times(mod):
    '''
    Find reverberation traveltimes from 1D model with taup
    '''
    model_in = mod.split('.tvel')[0]
    model = TauPyModel(model=model_in)
    phase_list = [
                  'ScSScS',
                  'sScSScS',
                  'ScSScSScS',
                  'sScSScSScS'
                  ]
    arrivals = model.get_travel_times(source_depth_in_km=160.,
                                   distance_in_degree=50.,
                                   phase_list = phase_list)
    time_array = np.zeros(5)
    time_array[0] = model_in.split('_')[1]

    for ii in range(0,len(arrivals)):
        time_array[1+ii]= np.round(arrivals[ii].time,1)

    return time_array

def diff_times(prem_dict,perturb_dict):
    '''
    Get traveltime delay between prem and perturbed model
    '''
    diff_dict = {}
    for key in prem_dict.keys():
        del_t = prem_dict[key] - perturb_dict[key]
        diff_dict[key] = del_t

    return diff_dict

tvel_list = glob('tvel/*.tvel')
prem_dict = reverb_times('prem')

for mod in tvel_list:
    print reverb_times(model)






