# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 14:52:47 2018

@author: Shanedrum
"""
import sys
sys.path.append("..")

import numpy as np
from astropy import constants as const
from astropy import units as u
from signals.ModSignal import *
from filters.RCHPF import *
from analysis.bin_consolidator import bin_consolidator
import h5py
from toolbox.axion_power import axion_power
from experiment.calc_sys_temp_offline import *



def convolve_two_arrays(array1, array2):
	"""Convolution based off the convolution theorem"""
	array1_fft = reciprocated_clone_hpf.DFT(array1)
	array2_fft = reciprocated_clone_hpf.DFT(array2)

	return reciprocated_clone_hpf.IDFT(array1_fft*array2_fft)

	
	
def MR_scan_analysis(scan_number, ch, res, modulation_type, axion_shape, notes):
    
	if not ch==2:
		ch=1
	file = h5py.File('run1a_data', 'r+')
	axion_scan = file["axion_log_run1a"][str(scan_number)]
	digitizer_scan = file["digitizer_log_run1a"][str(scan_number)]

	#Write exceptions here (reason not to include scans in Run1A)
	fstart = digitizer_scan.attrs["start_frequency"]
	fstop = digitizer_scan.attrs["stop_frequency"]
	res = digitizer_scan.attrs["frequency_resolution"]
	freqs = np.asarray([fstart + res*i for i in np.arange(start=0, stop=((fstop-fstart)/res))])

	binwidth = float(digitizer_scan.attrs["frequency_resolution"])*10**6
	h = const.h.value*u.J.to(u.eV,1) #plancks const eV*s
	k = const.k_B.value #Boltzmanns constant J/K
	Tsys = calc_sys_temp(digitizer_scan, axion_scan)["system temperature"] #calculate the system temperature
	kT = kTsys 
	BkT = binwidth*kT

	data = digitizer_scan[...]

	scan_length = len(freqs)
	middlefreqpos = np.floor(scan_length/2)
	middlefreq = digitizer_scan.freq[middlefreqpos-1]
	axion_coupling = 0.36
	dfszaxion = axion_power(ch, scan_number, None, axion_coupling, 0.45)
	axion_RMF = float(middlefreq*10**6 - (middlefreq*10**6)%binwidth) #calculate axion at center of scan
	timestamp = digitizer_scan["timestamp"]
	axion_mass = float(axion_RMF*h)
	startfreq = digitizer_scan.attrs["start_frequency"]*10**6
	wantTseries=None
					  
	signal_shape = modulatedsignal(modulation_type, timestamp, pCDM_w_baryons, axion_mass, resolution=binwidth)
	axblank = np.empty_like(signal_shape["signal"])
	DFSZshape = signal_shape*dfszaxion
	ax_freq = axion_mass/h

	#Remove Receiver response from scan
	filtered_data = reciprocated_clone_hpf(data, 3, ch)["filtereddata"]
	filtered_data_mean = np.mean(filtereddata)
	deltas = filtered_data - filtered_data_mean

	digitizer_scan = bin_consolidator(digitizer_scan, res)



	lorentzian_profile = axion_scan.lorentzian
	cc = 0.5
	cav_trans_mod = cc*lorentzian_profile
	axion_power_excess_watts = dfszaxion*cav_trans_mod

	#Genereate bin-wise scan stats assuming all power in single bin
	sigma = np.std(deltas)
	sigma_w = BkT*sigma
	power_deltas = BkT*deltas
	trans_power_deltas = cav_trans_mod*power_deltas
	variance_w = sigma_w**2
	nscans = lorentzian
	noise_power = sigma_w
	SNR = axion_power_excess_watts/sigma_w

	#Fit to axion signal

	#perform convolution based on convolution theorem

	trans_power_deltas_fft = reciprocated_clone_hpf.DFT(trans_power_deltas)
	DFSZshape_fft = reciprocated_clone_hpf.DFT(DFSZshape)
	dfszpow = dfszaxion
	signal_data_convolution = reciprocated_clone_hpf.IDFT(trans_power_deltas_fft*DFSZshape_fft)
	convolution_length = len(signal_data_convolution)

	axion_rmfs = []
	n_signal_width = len(axion_signal_shape.freqs)
	for i in np.arange(scan_length+2*n_signal_width)-1:
		axion_rmfs[i] = axion_RMF + binwidth*(i-middlefreqpos-n_signal_width)

	#find maximum likelihood

	conv_Aml_num = signal_data_convolution
	Aml_num = signal_data_convolution*(1/(2*sigma_w**2))
	lorentz_squared = cav_trans_mod**2
	model_squared = DFSZshape**2
	conv_Aml_den = convolve_two_arrays(lorentz_squared, model_squared)
	Aml_den = conv_Aml_den*(1/(2*sigma_w**2))
	model_excess_sqrd = Aml_den
	optimal_weight_sum = Aml_num

	#compute most likely axion power
	maximum_likelihood - Aml_num/Aml_den

	#Compute chi squared for maximum_likelihood sigma by explicit calculation
	A = maximum_likelihood
	array1 = [1]*len(A)
	A_factor = (2*A + 1)
	data_model_coeff = 2*Aml_num

	chi_sqrd_diff_one = array1*data_model_coeff
	chi_sqrd_diff_two = A_factor*Aml_den

	data_squared = power_deltas**2
	first_two_terms = A**2*Aml_den-2*A*Aml_num

	chi_squared_diff = chi_sqrd_diff_two - chi_sqrd_diff_one
	chi_squared_diff_red = chi_squared_diff/np.sqrt(len(data)-1)

	sigma_A = []
	for i in np.arange(len(chi_squared_diff_red)):
		sigma_A[i] = np.sqrt(2*chi_squared_diff[i])**(-1)

	#Fit significance
	axion_fit_significance = A/sigma_A

	#Sensitivity is calculated using significance of fit to axion of known power
	cl_coeff = 1.64485 #how many sigmas does 90% of data fall within
	sensitivity_power = sigma_A*cl_coeff
	sensitivity_coupling = np.sqrt(sensitivity_power)


	#consolidate statisitics
	nscans = np.append(axblank, np.append(nscans,axblank))
	SNR = np.append(axblank, np.append(SNR, axblank))
	noise_power = np.append(axblank, np.append(noise_power, axblank))
	weighted_deltas = np.append(axblank, np.append(weighted_deltas, axblank))
	power_deviation = np.append(axblank, np.append(power_deltas, axblank))


	results = {'scan_id':scan_number,
			   'nscans':nscans,
			   'sigma_w':sigma_w,
			   'optimal_weight_sum': optimal_weight_sum, #maximum likelihood numerator
			   'SNR':SNR,
			   'noise_power':noise_power,
			   'weighted_deltas':weighted_deltas,
			   'model_excess_sqrd':model_excess_sqrd, #maximum likelihood denominator
			   'axion fit':A,
			   'axion_fit_uncertainty':sigma_A,
			   'axion_fit_significance':axion_fit_significance,
			   'sensitivity_power':sensitivity_power,
			   'sensitivity_coupling':sensitivity_coupling,
			   'axion_frequencies':axion_rmfs,
			   'power deviation':power_deviation
			   }
	return results

        
    
    
    
    
    
    