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
from signals.modulation import *
from filters.RCHPF import *
from analysis.bin_consolidator import bin_consolidator
import h5py
from toolbox.axion_power import axion_power
from experiment.calc_sys_temp_offline import *
from toolbox.lorentzian import lorentzian



def convolve_two_arrays(array1, array2):
	"""Convolution based off the convolution theorem"""
	if len(array1)<len(array2):
		diff = len(array2)-len(array1)
		array1 = np.pad(array1, (int(np.floor(diff/2)), int(np.ceil(diff/2))), 'constant', constant_values=0)
	elif len(array2)<len(array1):
		diff = len(array1)-len(array2)
		array2 = np.pad(array2, (int(np.floor(diff/2)), int(np.ceil(diff/2))), 'constant', constant_values=0)
		
	array1_fft = DFT(array1)
	array2_fft = DFT(array2)

	return np.real(IDFT(array1_fft*array2_fft))

	
	
def MR_scan_analysis(scan,**params):
	"""
	Single scan analysis procedure.
	scanparams = dig_dataset, res, notes
	modparams = pec_vel, signal
	"""
	ch=1 #This is obsolete for 1A analysis, but was carried over from lua->python conversion.
	
	digitizer_scan = scan 
	scan_number = params["scan_number"]
	#Write exceptions here (reason not to include scans in Run1A)
	fstart = digitizer_scan.attrs["start_frequency"]
	fstop = digitizer_scan.attrs["stop_frequency"]
	res = digitizer_scan.attrs["frequency_resolution"]
	freqs = np.asarray([fstart + res*i for i in np.arange(start=0, stop=((fstop-fstart)/res))])

	binwidth = float(res)*10**6 # in Hz
	h = const.h.value*u.J.to(u.eV,1) #plancks const eV*s
	k = const.k_B.value #Boltzmanns constant J/K
	Tsys = params["Tsys"] #calculate the system temperature
	kT = k*Tsys 
	BkT = binwidth*kT

	data = digitizer_scan[...]

	modulation_type = params["pec_vel"]
	signal_shape = params["signal"]
	
	scan_length = len(freqs)
	middlefreqpos = int(np.floor(scan_length/2))
	middlefreq = freqs[(middlefreqpos-1)]
	nbins = params['nbins']
	axion_RMF = float(middlefreq*10**6 - (middlefreq*10**6)%binwidth) #calculate axion at center of scan
	timestamp = digitizer_scan.attrs["timestamp"]
	axion_mass = float(axion_RMF*h)
	startfreq = digitizer_scan.attrs["start_frequency"]*10**6
	wantTseries=None
	
	#Calculate average power deposited by axion
	dfszaxion = axion_power(params["axion_scan"],axion_RMF, nbins=nbins, res=res)
	
	#Calculate boosted signal
	mod_class = modulation()
	kwargs = {'resolution':binwidth}
	signal_shape = mod_class.modulatedsignal(modulation_type, timestamp, signal_shape, axion_mass, **kwargs)
	
	axblank = np.empty_like(signal_shape["signal"])
	DFSZshape = [i*dfszaxion for i in signal_shape["signal"]]

	#Remove Receiver response from scan
	filtered_data = reciprocated_clone_hpf(data, params["filter_params"][1])["filtereddata"]
	filtered_data_mean = np.mean(filtered_data)
	deltas = filtered_data - filtered_data_mean
	digitizer_scan = bin_consolidator(digitizer_scan, res)


	Q = eval(params["axion_scan"].attrs["Q"])
	res_freq = eval(params["axion_scan"].attrs["mode_frequency"])
	lorentzian_profile = lorentzian(Q, res_freq*10**6, nbins, binwidth)
	cc = 0.5
	cav_trans_mod = cc*lorentzian_profile
	axion_power_excess_watts = dfszaxion*cav_trans_mod

	#Genereate bin-wise scan stats assuming all power in single bin
	sigma = np.std(deltas)
	sigma_w = BkT*sigma
	power_deltas = BkT*deltas
	trans_power_deltas = cav_trans_mod*power_deltas
	variance_w = sigma_w**2
	nscans = lorentzian_profile
	noise_power = sigma_w
	SNR = axion_power_excess_watts/sigma_w

	#Fit to axion signal

	#perform convolution based on convolution theorem

	trans_power_deltas_fft = DFT(trans_power_deltas)
	DFSZshape_fft = DFT(DFSZshape)
	dfszpow = dfszaxion
	signal_data_convolution = convolve_two_arrays(trans_power_deltas, DFSZshape) 
	convolution_length = len(signal_data_convolution)

	axion_rmfs = []
	n_signal_width = len(DFSZshape)
	for i in np.arange(scan_length+2*n_signal_width)-1:
		axion_rmfs.append(axion_RMF + binwidth*(i-middlefreqpos-n_signal_width))

	#find maximum likelihood

	conv_Aml_num = signal_data_convolution
	Aml_num = signal_data_convolution*(1/(2*sigma_w**2))
	lorentz_squared = cav_trans_mod**2
	model_squared = [i**2 for i in DFSZshape] #DFSZshape**2
	conv_Aml_den = convolve_two_arrays(lorentz_squared, model_squared)
	Aml_den = conv_Aml_den*(1/(2*sigma_w**2))
	model_excess_sqrd = Aml_den
	optimal_weight_sum = Aml_num

	#compute most likely axion power
	maximum_likelihood = Aml_num/Aml_den
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
		sigma_A.append(np.sqrt(2*chi_squared_diff[i])**(-1))

	#Fit significance
	axion_fit_significance = A/sigma_A

	#Sensitivity is calculated using significance of fit to axion of known power
	cl_coeff = 1.64485 #how many sigmas does 90% of data fall within
	sensitivity_power = [i*cl_coeff for i in sigma_A] #sigma_A*cl_coeff
	sensitivity_coupling = [np.sqrt(i) for i in sensitivity_power] #np.sqrt(sensitivity_power)


	#consolidate statisitics
	nscans = np.append(axblank, np.append(nscans,axblank))
	SNR = np.append(axblank, np.append(SNR, axblank))
	noise_power = np.append(axblank, np.append(noise_power, axblank))
	power_deviation = np.append(axblank, np.append(power_deltas, axblank)) #weighted_deltas

	results = {'deltas':deltas,
				'scan_id':scan_number,
				'nscans':nscans,
				'sigma_w':sigma_w,
				'optimal_weight_sum': optimal_weight_sum, #maximum likelihood numerator
				'SNR':SNR,
				'noise_power':noise_power,
				'model_excess_sqrd':model_excess_sqrd, #maximum likelihood denominator
				'axion_fit':A,
				'axion_fit_uncertainty':sigma_A,
				'axion_fit_significance':axion_fit_significance,
				'sensitivity_power':sensitivity_power,
				'sensitivity_coupling':sensitivity_coupling,
				'axion_frequencies':axion_rmfs, #in Hz
				'power_deviation':power_deviation,
				'sigma':sigma
				}
	return results

        
    
    
    
    
    
    