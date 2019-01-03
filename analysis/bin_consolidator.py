# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 12:57:46 2018

@author: Shanedrum
Update 15-11-18: Updated script to handle h5py files
"""

import numpy as np
import h5py


def bin_consolidator(digitizer_dataset, new_resolution):
    """ 
    Consolidates power spectra into new resolution.
    Scan input is digitizer scan as pulled from receiver database, and new scan resolution in MHz
    bin counting starts at 0
	update 05-11-18: Configured for h5py file handling
    """
    scan = digitizer_dataset
    if not new_resolution:
        return "Error: new_resolution was not provided"
    elif int(new_resolution)<=0:
        return "Error: new_resolution must be positive"
    start_freq = scan["start_frequency"]
	stop_freq = scan["stop_frequency"]
	res = scan["frequency_resolution"]
	freqs = np.asarray([start_freq + res*i for i in np.arange(start=0, stop=((stop_freq-start_freq)/res))])
    #unpack scan
    freqs_ch1 = freqs
    freqs_ch2 = freqs #ch2 defaults to ch1
    power_spectrum_ch1 = scan[...]
    power_spectrum_ch2 = power_spectrum_ch1
    
    #divide old bins into new bins
    rebinned_freqs_ch1 = []
    rebinned_freqs_ch1.append(freqs_ch1[0] - freqs_ch1[0]%new_resolution)
    rebinned_freqs_ch2 = []
    rebinned_freqs_ch2.append(freqs_ch2[0] - freqs_ch2[0]%new_resolution)
    rebinned_power_spectrum_ch1 = np.empty_like(power_spectrum_ch1)
    rebinned_power_spectrum_ch2 = np.empty_like(power_spectrum_ch2)
    
    #ch1
    for i, val in enumerate(freqs_ch1):
        length = len(rebinned_freqs_ch1)
        length2 = len(freqs_ch1)
        adj_freq = val
        if i<length2-1:
            adj_freq = freqs_ch1[i+1]
        else:
            adj_freq = freqs_ch1[i-1]
            
        res = val-adj_freq
        if np.abs(val-rebinned_freqs_ch1[length])<= new_resolution/2.0:
            rebinned_power_spectrum_ch1[i] += power_spectrum_ch1[i]/res
        else:
            rebinned_freqs_ch1.append(rebinned_freqs_ch1[length]+new_resolution)
            rebinned_power_spectrum_ch1.append(power_spectrum_ch1[i]/res)
    rebinned_power_spectrum_ch1 = rebinned_power_spectrum_ch1/new_resolution
    
    #ch2
    for i, val in enumerate(freqs_ch2):
        length = len(rebinned_freqs_ch2)
        length2 = len(freqs_ch2)
        
        adj_freq = val
        if i<length2-1:
            adj_freq = freqs_ch2[i+1]
        else:
            adj_freq = freqs_ch2[i-1]
            
        res = val - adj_freq
        if np.abs(val - rebinned_freqs_ch2[length])<=new_resolution/2.0:
            rebinned_power_spectrum_ch2[i] += power_spectrum_ch2[i]/res
        else:
            rebinned_freqs_ch2.append(rebinned_freqs_ch2[length]+new_resolution)
            rebinned_power_spectrum_ch2.append(power_spectrum_ch2[i]/res)
    rebinned_power_spectrum_ch2 = rebinned_power_spectrum_ch2/new_resolution
    
    #deltas
    if scan.attrs["deltas"] and scan.attrs["frequencies_deltas"]:
        freqs_deltas = scan.attrs["frequencies_deltas"]
        deltas = scan.attrs["deltas"]
        rebinned_freqs_deltas = []
        rebinned_freqs_deltas[0] = freqs_deltas[0] - freqs_deltas[0]%new_resolution
        rebinned_deltas = []
        
        for i, val in enumerate(freqs_deltas):
            length = len(rebinned_freqs_deltas)
            length2 = len(freqs_deltas)
            adj_freq = val
            if i<length2-1:
                adj_freq = freqs_deltas[i+1]
            else:
                adj_freq = freqs_deltas[i-1]
        
            res = val - adj_freq
            if np.abs(val-rebinned_freqs_deltas[length])<= new_resolution/2.0:
                rebinned_deltas[i] += deltas[i]/res
            else:
                rebinned_freqs_deltas.append(rebinned_freqs_deltas[length] + new_resolution)
                rebinned_deltas.append(deltas[i]/res)
        rebinned_deltas = rebinned_deltas/new_resolution
    
    #update scan with rebinned quantities
    scan.attrs["start_frequency"] = rebinned_freqs_ch1[0]
	scan.attrs["stop_frequency"] = rebinned_freqs_ch1[-1]
	scan.attrs["frequency_resolution"] = new_resolution
    #scan.frequencies_channel_two = rebinned_freqs_ch2
    if scan.attrs["deltas"]:
        scan.attrs["frequencies_deltas"] = rebinned_freqs_deltas
    scan[...] = rebinned_power_spectrum_ch1
    #scan.power_spectrum_channel_two = rebinned_power_spectrum_ch2
    if scan.attrs["deltas"]:
        scan.attrs["deltas"] = rebinned_deltas
    
    return scan
            

    
    
    
    
    
    
    
    