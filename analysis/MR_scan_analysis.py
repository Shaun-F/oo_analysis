# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 14:52:47 2018

@author: Shaun Fell
"""
import sys
import os

import numpy
from astropy import constants as const
from astropy import units as U
import oo_analysis.filters.backsub_filters_lib as backsub_filters_lib
from oo_analysis.filters.backsub_filters_lib import *
from oo_analysis.toolbox.DFT import DFT, IDFT
from oo_analysis.analysis.bin_consolidator import bin_consolidator
import h5py
from oo_analysis.toolbox.axion_power import axion_power
from oo_analysis.experiment.calc_sys_temp_offline import *
from oo_analysis.toolbox.lorentzian import lorentzian
import time
import numba
from numba import jit
from scipy import fftpack
from scipy.signal.signaltools import _centered
from oo_analysis.analysis import synthetic_injection
from oo_analysis.analysis.synthetic_injection import axion_injector


def MR_scan_analysis(scan, **params):
	"""
	Single scan analysis procedure.
	Params must include: restype, notes, nbins, axion_frequencies_to_inject,
						pec_vel, signal_dataset, filter, filter_params, signal
	"""
	ch=1 #This is obsolete for 1A analysis, but was carried over from lua->python conversion.
	
	#declare variables to be used by single scan analysis
	try:
		digitizer_scan = scan 
		scan_number = params["scan_number"]
		#Write exceptions here (reason not to include scans in Run1A)
		
		fstart = float(digitizer_scan.attrs["start_frequency"])
		fstop = float(digitizer_scan.attrs["stop_frequency"])
		res = float(digitizer_scan.attrs["frequency_resolution"])
		int_time = float(digitizer_scan.attrs['integration_time'])
		freqs = numpy.asarray([fstart + res*i for i in numpy.arange(start=0, stop=((fstop-fstart)/res))])

		binwidth = float(res)*10**6 # in Hz
		
		
				
		modulation_type = params["pec_vel"]
		signal_shape = params["signal_dataset"][scan_number]['signal']
		signal_width_Hz = params["signal_dataset"][scan_number]['signal_width']
		data = digitizer_scan[...]
		
		h = const.h.to(U.eV*U.s).value #plancks const eV*s
		k = const.k_B.to(U.J/U.Kelvin).value #Boltzmanns constant J/K
		Tsys = params["Tsys"] #calculate the system temperature
		power_johnson = signal_width_Hz*k*Tsys #johnson noise power
		
		scan_length = len(freqs)
		middlefreqpos = int(numpy.floor(scan_length/2))
		middlefreq = freqs[(middlefreqpos-1)]
		nbins = params['nbins']
		axion_RMF = float(middlefreq*10**6 - (middlefreq*10**6)%binwidth) #calculate axion at center of scan (Hz)
		timestamp = digitizer_scan.attrs["timestamp"]
		axion_mass = float(axion_RMF*h)
		startfreq = float(digitizer_scan.attrs["start_frequency"])*10**6 #Hz
		wantTseries=None
		
		
		
		
		#Calculate average power deposited by axion
		dfszaxion = axion_power(digitizer_scan,axion_RMF)
		
	except (SyntaxError, KeyError) as error:
		print("MR_scan_analysis failed at scan {0} with error: \n {1}".format(scan_number, error))
		raise
	
	#begin signal scan analysis
	try:		
		axblank = numpy.empty_like(signal_shape)
		DFSZshape = signal_shape*dfszaxion #Units of watts

		#Remove Receiver response from scan

		filter_function = getattr(backsub_filters_lib, params['filter'])
		data[0] = np.mean([data[1], data[2], data[3]]) # Some spectra have odd behavior at beginning of scan, i.e. a single downward spike at the beginning position. I just set a default value
		filtered_dict = filter_function(data, **{**params['filter_params'], **params})
		filtered_data = filtered_dict['filtereddata']
	
		#Inject artificial axion
		
		params['fstart'] = fstart
		params['fstop'] = fstop
		params['fres'] = res
		
		#Inject software-generated axion signals
		synthetic_injector = axion_injector(digitizer_scan, filter_shape=filtered_dict['background'].real, **params)
		data, SIA_meta = synthetic_injector.inject_signal()
		filtered_dict = filter_function(data, **params['filter_params'])
		filtered_data = filtered_dict['filtereddata']

		
		filtered_data_mean = numpy.mean(filtered_data)
		deltas = np.asarray((filtered_data - filtered_data_mean))
		
		#digitizer_scan = bin_consolidator(digitizer_scan, res)
		
		
		
		Q = digitizer_scan.attrs["Q"]
		res_freq = digitizer_scan.attrs["mode_frequency"] #MHz
		lorentzian_profile = lorentzian(Q, res_freq*10**6, fstart*10**6, fstop*10**6, binwidth)
		cc = 0.5
		cav_trans_mod = cc*lorentzian_profile
		axion_power_excess_watts = convolve_two_arrays(DFSZshape, cav_trans_mod)

		
		#Genereate bin-wise scan stats assuming all power in single bin
		sigma = numpy.std(deltas)
		sigma_w = power_johnson*(signal_width_Hz*int_time)**(-1/2)
		power_deltas = power_johnson*deltas
		nscans = lorentzian_profile
		SNR = (axion_power_excess_watts/sigma_w)
		
		"""
		#Testing new SNR calculation
		WIENER = np.convolve(deltas, DFSZshape)
		SIGMA = np.sqrt(sigma**2 * sum(DFSZshape**2))
		SNR = WIENER/SIGMA
		"""
		
		chi_squared_results = chi_squared(power_deltas, DFSZshape, lorentzian_profile,  sigma_w, cc=0.5)
		
		sensitivity_power = chi_squared_results['power_sensitivity']
		sensitivity_coupling = chi_squared_results['coupling_sensitivity']
		maximum_likelihood = chi_squared_results['maximum_likelihood']
		axion_fit_significance = chi_squared_results['axion_fit_significance']
		axion_fit_uncertainty = chi_squared_results['axion_fit_uncertainty']
		optimal_weight_sum = chi_squared_results['chi_squared_term_two']
		model_excess_sqrd = chi_squared_results['chi_squared_term_three']
				
		#Candidate flagging
		
		
		axion_rmfs = []
		n_signal_width = len(DFSZshape)
		for i in numpy.arange(scan_length+2*n_signal_width)-1:
			axion_rmfs.append(axion_RMF + binwidth*(i-middlefreqpos-n_signal_width))
		

		#consolidate statisitics
		nscans = numpy.pad(nscans, len(axblank), 'constant', constant_values = 0)
		SNR = numpy.pad(SNR, len(axblank), 'constant', constant_values = 0)
		
		
		
	except (KeyError, ValueError, IndexError) as error:
		print("\n\nError with scan {0} in single scan analysis script.".format(scan_number))
		open(os.getcwd() + "/oo_analysis/meta/error_log", 'a+').write("\n\n"+ str(error))
		raise
	
	results = {'deltas':deltas,
				'scan_id':scan_number,
				'nscans':nscans,
				'sigma_w':sigma_w,
				'optimal_weight_sum': optimal_weight_sum, #maximum likelihood numerator
				'SNR':SNR,
				'noise_power':sigma_w,
				'model_excess_sqrd':model_excess_sqrd, #maximum likelihood denominator
				'axion_fit':maximum_likelihood,
				'axion_fit_uncertainty':axion_fit_uncertainty,
				'axion_fit_significance':axion_fit_significance,
				'sensitivity_power':sensitivity_power,
				'sensitivity_coupling':sensitivity_coupling,
				'power_deviation':power_deltas.real,
				'sigma':sigma,
				'start_frequency': freqs[0]*10**6,
				'middle_frequency':middlefreq*10**6,
				'axion_frequencies':axion_rmfs, #in Hz
				'software_injected_axion_properties':SIA_meta
				}
	return results

def convolve_two_arrays(array1, array2):
	"""Convolution based off the convolution theorem"""
	if len(array1)<len(array2):
		diff = len(array2)-len(array1)
		array1 = numpy.pad(array1, (int(numpy.floor(diff/2)), int(numpy.ceil(diff/2))), 'constant', constant_values=0)
	elif len(array2)<len(array1):
		diff = len(array1)-len(array2)
		array2 = numpy.pad(array2, (int(numpy.floor(diff/2)), int(numpy.ceil(diff/2))), 'constant', constant_values=0)
		
	#The next 8 lines follow closely with scipy.signal.fftconvolve
	shape = np.maximum(np.array(array1.shape), np.array(array2.shape))
	s1 = np.array(array1.shape)
	s2 = np.array(array2.shape)
	shape = s1+s2-1
	fshape = [fftpack.helper.next_fast_len(d) for d in shape]
	fslice = tuple([slice(sz) for sz in shape])
	array1_fft = DFT(array1, fshape)
	array2_fft = DFT(array2, fshape)
	
	convolved_arr = numpy.asarray(numpy.real(IDFT(array1_fft*array2_fft)))[fslice]
	ret = _centered(convolved_arr, s1)
	return ret

def chi_squared(deltas, axion_signal, cavity_lorentzian, noise_disp, convolve=True, **kwargs):
	"""
	Chi squared as calculated in Improving Axion Signal Models Through N-Body 
	Simulations by Erik Lentz 2017. Specifically, section 5.5. Each term corresponds to the terms of equation 5.23 in thesis.
	NOTE: Axion_signal should be in units of watts
	"""
    
	if 'cc' not in kwargs:
		cc=0.5  #coupling of photons to cavity
	else:
		cc = kwargs['cc']
	if 'cl_coeff' not in kwargs:
		cl_coeff = 1.644853 # For 90% Confidence interval, set z = 1.644853  For 95% C.I., set z = 1.959964
	else:
		cl_coef = kwargs['cl_coeff']
	pad_len = len(axion_signal)
	
	cavity_transmission = cc*cavity_lorentzian #tranmission function for cavity
	transmitted_power_deltas = numpy.pad(cavity_transmission*deltas, pad_len, 'constant', constant_values = 0) #power deltas are dimensionless deltas times Bandwidth*kT
	cavity_transmission = numpy.pad(cavity_transmission, pad_len, 'constant', constant_values=0)
	
	
	
	chi_squared_term_one = numpy.sum(deltas**2/(2*noise_disp**2))
	chi_squared_term_two = convolve_two_arrays(transmitted_power_deltas, axion_signal)/(2*noise_disp**2)
	chi_squared_term_three = convolve_two_arrays(cavity_transmission**2, axion_signal**2)/(2*noise_disp**2)
	
	chi_squared = lambda A: chi_squared_term_one - 2*A*chi_squared_term_two + A**2 * chi_squared_term_three
	
	maximum_likelihood = chi_squared_term_two/chi_squared_term_three #equ 5.26 in Lentz Thesis
	
	chi_squared_small_difference = chi_squared(maximum_likelihood+1) - chi_squared(maximum_likelihood) #equation 5.28 in Lentz Thesis
	
	maximum_likelihood_uncertainty = (2*chi_squared_small_difference)**(-1/2) #equ 5.27 in Lentz Thesis
	
	maximum_likelihood_uncertainty = new_padding(maximum_likelihood_uncertainty, (pad_len,pad_len), pad_val = np.inf) 
	
	axion_fit_significance = maximum_likelihood/maximum_likelihood_uncertainty#equ 5.34 in Lentz Thesis
	
	#Sensitivity is calculated using significance of fit to axion of known power
	power_sensitivity = new_padding(maximum_likelihood_uncertainty*cl_coeff, (pad_len,pad_len), pad_val = np.inf) #Pow. Sens goes like 1/axion power
	coupling_sensitivity = numpy.sqrt(power_sensitivity)

	"""
	print(numpy.nanstd(axion_fit_significance))
	import matplotlib.pyplot as plt
	fig, ax = plt.subplots(4,2)
	ax[0,0].title("power Deltas")
	ax[0,0].plot(deltas)
	ax[0,1].title("power deltas histogram")
	ax[0,1].hist(deltas, bins=25, align='mid', histtype='step')
	ax[0,2].title("axion fit sign")
	ax[0,2].plot(axion_fit_significance)
	ax[0,3].title("maximum likelihood")
	ax[0,3].plot(maximum_likelihood)
	ax[0,4].title("maximum likelihood dispersion")
	ax[0,4].plot(maximum_likelihood_uncertainty)
	ax[1,0].title("chi_square small diff")
	ax[1,0].plot(chi_squared_small_difference)
	ax[1,1].title("Coupling Sensitivity")
	ax[1,1].plot(coupling_sensitivity)
	plt.show()
	"""
	
	results = {'chi_squared_term_one': chi_squared_term_one,
				'chi_squared_term_two': chi_squared_term_two,
				'chi_squared_term_three': chi_squared_term_three,
				'maximum_likelihood': maximum_likelihood,
				'axion_fit_significance': axion_fit_significance,
				'power_sensitivity': power_sensitivity,
				'coupling_sensitivity': coupling_sensitivity,
				'axion_fit_uncertainty': maximum_likelihood_uncertainty
				}
	
	return results
	
	
def flip_arr(arr):
	mid = int(numpy.ceil(len(arr)/2))
	new_arr = numpy.empty_like(arr)
	
	if len(arr)%2==0: #If length of array is even
		new_arr[:mid] = arr[mid:]
		new_arr[mid:] = arr[:mid]
	else: #if length of array is odd
		new_arr[:mid] = arr[mid-1:]
		new_arr[mid:] = arr[:mid-1]
	return numpy.asarray(new_arr)
	

def remove_padding(arr, pad_len): #pad_len is tuple-like. zeroth index is length of left padding. first index is length of right padding
	new_arr = numpy.empty_like(arr)
	arr_len = len(arr)
	lpad, rpad = pad_len[0], pad_len[1]

	new_arr = arr[lpad:arr_len-rpad]
	return new_arr

def new_padding(arr, pad_len, pad_val = 0): #pad_len is tuple-like. zeroth index is length of left padding. first index is length of right padding
	new_arr = remove_padding(arr, pad_len)
	lpad, rpad = pad_len[0], pad_len[1]
	new_arr = numpy.pad(new_arr, (lpad, rpad), 'constant', constant_values=pad_val)
	return new_arr
