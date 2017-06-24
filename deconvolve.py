#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : deconvolve.py
Purpose : deconvolve windows
Creation Date : 15-06-2017
Last Modified : Fri 23 Jun 2017 07:30:28 PM EDT
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
from scipy.linalg import toeplitz
from scipy.signal import cosine
from scipy.signal import gaussian
from scipy.signal import tukey
from scipy.signal import boxcar


def main():
    #st = read_stream('/home/samhaug/work1/ScS_reverb_data/obspyDMT/mag_7.0_8.0/20100812_115415.a/processed/')
    #st = obspy.read('/home/samhaug/work1/ScS_reverb_data/013016E/FJ.pk')
    #st = obspy.read('/home/samhaug/work1/ScS_reverb_sims/mineos/PREM_RUNS/prem_013016_japan/st_T.pk')
    #st = obspy.read('/home/samhaug/work1/ScS_reverb_data/obspyDMT/mag_7.0_8.0/20151124_224538.a/processed/st_clean_unfiltered.pk')
    st = obspy.read('/home/samhaug/work1/ScS_reverb_data/obspyDMT/mag_7.0_8.0/20100812_115415.a/processed/st_T_unfiltered.pk')
    st.interpolate(10)
    #st.differentiate()
    st.filter('bandpass',freqmin=1./100,freqmax=1./10,zerophase=True)
    st.normalize()
    seispy.plot.simple_section(st)
    st = seispy.filter.range_filter(st,(35,44))
    #st.differentiate()
    #st.differentiate().normalize()
    print st[0].stats.sac['evdp']

    a = []
    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    for idx,tr in enumerate(st):
        tr.stats.sac['o'] = 0
        f,mask = mask_trace(tr,['sScSScS'])
        ax2.plot(f.data+idx,color='k',alpha=0.3)
        ax2.plot(mask+idx,color='g')
        rf = water_level(f,mask,alpha=0.05,plot=False)
        rf = np.roll(rf,int(st[0].stats.sampling_rate*25)-np.argmax(rf))
        a.append(rf/rf.max())
        ax1.plot(rf/rf.max()+idx)
    plt.show()
    t = np.linspace(-50,500,f.stats.npts)
    plt.plot(t,np.mean(a,axis=0))
    plt.show()

def read_stream(dirname):
    st = obspy.read(dirname+'/st_T_clean.pk')
    #st.filter('highpass',freq=1./80,zerophase=True)
    return st

def mask_trace(tr,phase):
    sr = tr.stats.sampling_rate
    f = seispy.data.phase_window(tr,phase=phase,window=(-50,500))
    mask = np.zeros(len(f.data))
    mask[int(0*sr):int(100*sr)] = boxcar(int(100*sr))
    mask = f.data*mask
    return f,mask

def water_level(a,b,alpha=0.1,plot=False):
    t = np.linspace(0,a.stats.npts*a.stats.sampling_rate,a.stats.npts)

    #Convert to frequency domain-------------------------------
    a_omega = np.fft.fft(a.data)
    b_omega = np.fft.fft(b)

    #Perform division------------------------------------------
    F_omega = b_omega*b_omega.conjugate()
    Phi_ss  = np.maximum(F_omega,alpha*(np.amax(F_omega)))
    try:
        H_omega = ((a_omega*b_omega.conjugate())/Phi_ss)
    except RuntimeWarning:
        return np.zeros(len(a.data))

    #Convert back to time domain-------------------------------
    #rf = np.zeros(len(H_omega))
    rf = np.fft.ifft(H_omega)

    #Plots-----------------------------------------------------
    if plot==True:
        fig,ax = plt.subplots(2,1)
        ax[0].plot(t,a)
        ax[0].plot(t,b)
        ax[1].plot(t,rf)
        plt.show()

    return np.real(rf)

main()



