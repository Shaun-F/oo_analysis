# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 14:52:47 2018

@author: Shaun Fell
"""
import sys
sys.path.append("..")

import numpy
from astropy import constants as const
from astropy import units as u
from signals.modulation import *
from filters.backsub_filters_lib import RCHPF as reciprocated_clone_hpf
from toolbox.DFT import DFT, IDFT
from analysis.bin_consolidator import bin_consolidator
import h5py
from toolbox.axion_power import axion_power
from experiment.calc_sys_temp_offline import *
from toolbox.lorentzian import lorentzian
import time
import numba; from numba import jit
from scipy import fftpack
from scipy.signal.signaltools import _centered


def MR_scan_analysis(scan, **params):
	"""
	Single scan analysis procedure.
	scanparams = dig_dataset, res, notes
	modparams = pec_vel, signal
	"""
	ch=1 #This is obsolete for 1A analysis, but was carried over from lua->python conversion.
	
	#declare variables to be used by single scan analysis
	try:
		digitizer_scan = scan 
		scan_number = params["scan_number"]
		#Write exceptions here (reason not to include scans in Run1A)
		
		constants_start = time.time()
		fstart = float(digitizer_scan.attrs["start_frequency"])
		fstop = float(digitizer_scan.attrs["stop_frequency"])
		res = float(digitizer_scan.attrs["frequency_resolution"])
		int_time = float(digitizer_scan.attrs['integration_time'])
		freqs = numpy.asarray([fstart + res*i for i in numpy.arange(start=0, stop=((fstop-fstart)/res))])

		binwidth = float(res)*10**6 # in Hz
		
		h = const.h.to(U.eV*u.s).value #plancks const eV*s
		k = const.k_B.to(U.J/U.Kelvin).value #Boltzmanns constant J/K
		Tsys = params["Tsys"] #calculate the system temperature
		power_johnson = binwidth*k*Tsys
				
		modulation_type = params["pec_vel"]
		signal_shape = params["signal_dataset"][scan_number]['signal']
		
		data = digitizer_scan[...]
		
		scan_length = len(freqs)
		middlefreqpos = int(numpy.floor(scan_length/2))
		middlefreq = freqs[(middlefreqpos-1)]
		nbins = params['nbins']
		axion_RMF = float(middlefreq*10**6 - (middlefreq*10**6)%binwidth) #calculate axion at center of scan (Hz)
		timestamp = digitizer_scan.attrs["timestamp"]
		axion_mass = float(axion_RMF*h)
		startfreq = float(digitizer_scan.attrs["start_frequency"])*10**6 #Hz
		wantTseries=None
		
		constants_stop = time.time()
		
		
		
		#Calculate average power deposited by axion
		axion_power_start = time.time()
		dfszaxion = axion_power(params["axion_scan"],axion_RMF, nbins=nbins, res=res)
		axion_power_stop = time.time()
		#print("dfszaxion", dfszaxion)
	except SyntaxError as error:
		print("MR_scan_analysis failed at scan {0} with error: \n {1}".format(scan_number, error))
		raise
	
	#begin signal scan analysis
	try:
		#Calculate boosted signal
		modulation_start = time.time()
		mod_class = modulation()
		kwargs = {'resolution':binwidth}
		modulation_stop = time.time()
		
		axblank = numpy.empty_like(signal_shape)
		DFSZshape = signal_shape*dfszaxion

		#Remove Receiver response from scan
		BS_start = time.time()
		try:
			copies=params["filter_params"][1]
			window = params['filter_params'][0]
			filter, errors = reciprocated_clone_hpf(data, window, copies, False, **params['submeta'])
			filtered_data = filter['filtereddata']
			submeta = filter['meta']
			#print("copies", copies, "\n window", window, "\n filtered_data", filtered_data)
			if errors['maxed_filter_size']:
				with open('BS_errors.txt', 'a+') as f:
					f.write("Background subtraction filter size maxed out with scan {0} \n".format(scan_number))
		except TypeError as error:
			raise
			#raise Exception("Error: Invalid data type for scan {0}. Type {1} not accepted".format(scan_number, type(data)))
		BS_stop = time.time()
		
		consolidation_start = time.time()
		filtered_data_mean = numpy.mean(filtered_data)
		deltas = np.asarray((filtered_data - filtered_data_mean))
		digitizer_scan = bin_consolidator(digitizer_scan, res)
		
		
		consolidation_stop = time.time()
		
		
		cavity_lorentz_start = time.time()
		Q = params["axion_scan"].attrs["Q"]
		res_freq = params["axion_scan"].attrs["mode_frequency"] #MHz
		lorentzian_profile = lorentzian(Q, res_freq*10**6, fstart*10**6, fstop*10**6, binwidth)
		cc = 0.5
		cav_trans_mod = cc*lorentzian_profile
		axion_power_excess_watts = convolve_two_arrays(DFSZshape, cav_trans_mod)
		cavity_lorentz_stop = time.time()

		
		#Genereate bin-wise scan stats assuming all power in single bin
		bin_stats_start = time.time()
		sigma = numpy.std(deltas)
		sigma_w = power_johnson*(binwidth*int_time)**(-1/2)
		power_deltas = power_johnson*deltas
		nscans = lorentzian_profile
		SNR = (axion_power_excess_watts/sigma_w)
		bin_stats_stop = time.time()
		
		
		chi_squared_start = time.time()
		chi_squared_results = chi_squared(power_deltas, DFSZshape, lorentzian_profile, sigma_w, cc=0.5)
		chi_squared_stop = time.time()
		
		sensitivity_power = chi_squared_results['power_sensitivity']
		sensitivity_coupling = chi_squared_results['coupling_sensitivity']
		maximum_likelihood = chi_squared_results['maximum_likelihood']
		axion_fit_significance = chi_squared_results['axion_fit_significance']
		axion_fit_uncertainty = chi_squared_results['axion_fit_uncertainty']
		optimal_weight_sum = chi_squared_results['chi_squared_term_two']
		model_excess_sqrd = chi_squared_results['chi_squared_term_three']
		#Fit to axion signal

		
		axion_rmfs_start = time.time()
		axion_rmfs = []
		n_signal_width = len(DFSZshape)
		for i in numpy.arange(scan_length+2*n_signal_width)-1:
			axion_rmfs.append(axion_RMF + binwidth*(i-middlefreqpos-n_signal_width))
		axion_rmfs_stop = time.time()
		

		#consolidate statisitics
		nscans = numpy.pad(nscans, len(axblank), 'constant', constant_values = 0)
		SNR = numpy.pad(SNR, len(axblank), 'constant', constant_values = 0)
		
		
		
	except (KeyError, ValueError, IndexError) as error:
		print("\n\nError with scan {0} in single scan analysis script.".format(scan_number))
		open('../meta/error_log', 'a+').write(str(time.time())+ "\n\n"+ str(error))
		raise
		
		
	if submeta['timeit']:
		submeta = params['submeta']
		submeta['constants'].append(constants_stop - constants_start)
		submeta['axion_power'].append(axion_power_stop - axion_power_start)
		submeta['modulation'].append(modulation_stop - modulation_start)
		submeta['BS'].append(BS_stop - BS_start)
		submeta['convolutions'].append(chi_squared_stop - chi_squared_start)
		submeta['consolidation'].append(consolidation_stop - consolidation_start)
		submeta['cavity_lorentz'].append(cavity_lorentz_stop - cavity_lorentz_start)
		submeta['bin_stats'].append(bin_stats_stop - bin_stats_start)
		submeta['axion_rmfs'].append(axion_rmfs_stop - axion_rmfs_start)
	
	
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
				'axion_frequencies':axion_rmfs #in Hz
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

def chi_squared(power_deltas, axion_signal, cavity_lorentzian, noise_power, convolve=True, **kwargs):
	"""
	Chi squared as calculated in Improving Axion Signal Models Through N-Body 
	Simulations by Erik Lentz 2017. Specifically, section 5.5. Each term corresponds to the terms of equation 5.23 in thesis.
	"""
    
	if 'cc' not in kwargs:
		cc=0.5  #coupling of photons to cavity
	else:
		cc = kwargs['cc']
	if 'cl_coeff' not in kwargs:
		cl_coeff = 1.64485 #90% confidence interval
	else:
		cl_coef = kwargs['cl_coeff']
	pad_len = len(axion_signal)
	
	
	
	cavity_transmission = cc*cavity_lorentzian #tranmission function for cavity
	transmitted_power_deltas = numpy.pad(cavity_transmission*power_deltas, pad_len, 'constant', constant_values = 0) #power deltas are dimensionless deltas times Bandwidth*kT
	cavity_transmission = numpy.pad(cavity_transmission, pad_len, 'constant', constant_values=0)
	
	chi_squared_term_one = (1/2)*numpy.sum(power_deltas**2/noise_power)
	chi_squared_term_two = convolve_two_arrays(transmitted_power_deltas, axion_signal)/(2*noise_power**2) #also having issues
	chi_squared_term_three = convolve_two_arrays(cavity_transmission**2, axion_signal**2)/(2*noise_power**2) #this is causing problems
	
	maximum_likelihood = chi_squared_term_two/chi_squared_term_three #equ 5.26 in Lentz Thesis
	
	chi_squared_small_difference = chi_squared_term_three*(2*maximum_likelihood+1) - 2*chi_squared_term_two #equation 5.28 in Lentz Thesis
	maximum_likelihood_dispersion = new_padding(1/numpy.sqrt(2*chi_squared_small_difference), (pad_len,pad_len), pad_val = np.inf) #equ 5.27 in Lentz Thesis
	
	axion_fit_significance = maximum_likelihood/maximum_likelihood_dispersion#equ 5.34 in Lentz Thesis
	
	#Sensitivity is calculated using significance of fit to axion of known power
	power_sensitivity = new_padding(maximum_likelihood_dispersion*cl_coeff, (pad_len,pad_len), pad_val = np.inf)
	coupling_sensitivity = numpy.sqrt(power_sensitivity)

	
	results = {'chi_squared_term_one': chi_squared_term_one,
				'chi_squared_term_two': chi_squared_term_two,
				'chi_squared_term_three': chi_squared_term_three,
				'maximum_likelihood': maximum_likelihood,
				'axion_fit_significance': axion_fit_significance,
				'power_sensitivity': power_sensitivity,
				'coupling_sensitivity': coupling_sensitivity,
				'axion_fit_uncertainty': maximum_likelihood_dispersion
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
	
	
	
	
	
	
	
	
	
    