# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 14:52:47 2018

@author: Shanedrum
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



def convolve_two_arrays(array1, array2):
	"""Convolution based off the convolution theorem"""
	if len(array1)<len(array2):
		diff = len(array2)-len(array1)
		array1 = numpy.pad(array1, (int(numpy.floor(diff/2)), int(numpy.ceil(diff/2))), 'constant', constant_values=0)
	elif len(array2)<len(array1):
		diff = len(array1)-len(array2)
		array2 = numpy.pad(array2, (int(numpy.floor(diff/2)), int(numpy.ceil(diff/2))), 'constant', constant_values=0)
		
	array1_fft = DFT(array1)
	array2_fft = DFT(array2)

	return numpy.real(IDFT(array1_fft*array2_fft))

	
	
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
		
		h = const.h.value*u.J.to(u.eV,1) #plancks const eV*s
		k = const.k_B.value #Boltzmanns constant J/K
		Tsys = params["Tsys"] #calculate the system temperature
		kT = k*Tsys 
		BkT = binwidth*kT
		constants_stop = time.time()
		
		data = digitizer_scan[...]

		
		modulation_type = params["pec_vel"]
		signal_shape = params["signal_dataset"][scan_number]['signal']
		
		scan_length = len(freqs)
		middlefreqpos = int(numpy.floor(scan_length/2))
		middlefreq = freqs[(middlefreqpos-1)]
		nbins = params['nbins']
		axion_RMF = float(middlefreq*10**6 - (middlefreq*10**6)%binwidth) #calculate axion at center of scan
		timestamp = digitizer_scan.attrs["timestamp"]
		axion_mass = float(axion_RMF*h)
		startfreq = float(digitizer_scan.attrs["start_frequency"])*10**6
		wantTseries=None
		
		#Calculate average power deposited by axion
		axion_power_start = time.time()
		dfszaxion = axion_power(params["axion_scan"],axion_RMF, nbins=nbins, res=res)
		axion_power_stop = time.time()
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
		DFSZshape = [i*dfszaxion for i in signal_shape]

		#Remove Receiver response from scan
		BS_start = time.time()
		try:
			copies=params["filter_params"][1]
			window = params['filter_params'][0]
			filter, errors = reciprocated_clone_hpf(data, window, copies, False, **params['submeta'])
			filtered_data = filter['filtereddata']
			submeta = filter['meta']
			if errors['maxed_filter_size']:
				with open('BS_errors.txt', 'a+') as f:
					f.write("Background subtraction filter size maxed out with scan {0} \n".format(scan_number))
		except TypeError as error:
			raise
			#raise Exception("Error: Invalid data type for scan {0}. Type {1} not accepted".format(scan_number, type(data)))
		BS_stop = time.time()
		
		consolidation_start = time.time()
		filtered_data_mean = numpy.mean(filtered_data)
		deltas = (filtered_data - filtered_data_mean)
		digitizer_scan = bin_consolidator(digitizer_scan, res)
		consolidation_stop = time.time()
		
		
		cavity_lorentz_start = time.time()
		Q = params["axion_scan"].attrs["Q"]
		res_freq = params["axion_scan"].attrs["mode_frequency"]
		lorentzian_profile = lorentzian(Q, res_freq*10**6, nbins, binwidth)
		cc = 0.5
		cav_trans_mod = cc*lorentzian_profile
		axion_power_excess_watts = dfszaxion*cav_trans_mod
		cavity_lorentz_stop = time.time()
		
		
		#Genereate bin-wise scan stats assuming all power in single bin
		bin_stats_start = time.time()
		sigma = numpy.std(deltas)
		sigma_w = BkT*sigma
		power_deltas = BkT*deltas
		trans_power_deltas = cav_trans_mod*power_deltas
		variance_w = sigma_w**2
		nscans = lorentzian_profile
		noise_power = sigma_w
		SNR = (axion_power_excess_watts/sigma_w)*(binwidth*int_time)**(1/2)
		bin_stats_stop = time.time()
		
		#Fit to axion signal

		#perform convolution based on convolution theorem
		
		convolutions_start = time.time()
		trans_power_deltas_fft = DFT(trans_power_deltas)
		DFSZshape_fft = DFT(DFSZshape)
		dfszpow = dfszaxion
		trans_power_deltas = numpy.pad(trans_power_deltas, len(axblank), 'constant', constant_values = 0)
		signal_data_convolution = convolve_two_arrays(trans_power_deltas, DFSZshape) 
		convolution_length = len(signal_data_convolution)
		convolutions_stop = time.time()
		
		
		axion_rmfs_start = time.time()
		axion_rmfs = []
		n_signal_width = len(DFSZshape)
		for i in numpy.arange(scan_length+2*n_signal_width)-1:
			axion_rmfs.append(axion_RMF + binwidth*(i-middlefreqpos-n_signal_width))
		axion_rmfs_stop = time.time()

		
		#find maximum likelihood

		max_likelihood_arith_start = time.time()
		conv_Aml_num = signal_data_convolution
		Aml_num = signal_data_convolution*(1/(2*sigma_w**2))
		lorentz_squared = cav_trans_mod**2
		model_squared = [i**2 for i in DFSZshape] #DFSZshape**2
		lorentz_squared = numpy.pad(lorentz_squared, len(axblank), 'constant', constant_values=0)
		conv_Aml_den = convolve_two_arrays(lorentz_squared, model_squared)
		Aml_den = conv_Aml_den*(1/(2*sigma_w**2))
		model_excess_sqrd = Aml_den
		optimal_weight_sum = Aml_num
		max_likelihood_arith_stop = time.time()
		#compute most likely axion power
		max_likelihood_start = time.time()
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
		chi_squared_diff_red = chi_squared_diff/numpy.sqrt(len(data)-1)

		sigma_A = []
		for i in numpy.arange(len(chi_squared_diff_red)):
			sigma_A.append(numpy.sqrt(2*chi_squared_diff[i])**(-1))
		max_likelihood_stop = time.time()
		
		sig_sens_start = time.time()
		#Fit significance
		axion_fit_significance = A/sigma_A

		
		#Sensitivity is calculated using significance of fit to axion of known power
		cl_coeff = 1.64485 #how many sigmas does 90% of data fall within
		sensitivity_power = [i*cl_coeff for i in sigma_A] #sigma_A*cl_coeff
		sensitivity_coupling = [numpy.sqrt(i) for i in sensitivity_power] #np.sqrt(sensitivity_power)


		#consolidate statisitics
		nscans = numpy.pad(nscans, len(axblank), 'constant', constant_values = 0)
		#deltas = numpy.pad(deltas, len(axblank), 'constant', constant_values=0)
		SNR = numpy.pad(SNR, len(axblank), 'constant', constant_values = numpy.inf)
		#power_deviation = numpy.pad(power_deltas, len(axblank), 'constant', constant_values = 0) #weighted_deltas
		sig_sens_stop = time.time()
	except (KeyError, ValueError, IndexError) as error:
		print("\n\nError with scan {0}".format(scan_number))
		open('../../meta/error_log', 'a+').write(str(time.time())+ "\n\n"+ str(error))
		raise
		
		
	if submeta['timeit']:
		submeta = params['submeta']
		submeta['constants'].append(constants_stop - constants_start)
		submeta['axion_power'].append(axion_power_stop - axion_power_start)
		submeta['modulation'].append(modulation_stop - modulation_start)
		submeta['BS'].append(BS_stop - BS_start)
		submeta['consolidation'].append(consolidation_stop - consolidation_start)
		submeta['cavity_lorentz'].append(cavity_lorentz_stop - cavity_lorentz_start)
		submeta['bin_stats'].append(bin_stats_stop - bin_stats_start)
		submeta['convolutions'].append(convolutions_stop - convolutions_start)
		submeta['axion_rmfs'].append(axion_rmfs_stop - axion_rmfs_start)
		submeta['max_likelihood_arith'].append(max_likelihood_arith_stop - max_likelihood_arith_start)
		submeta['max_likelihood'].append(max_likelihood_stop - max_likelihood_start)
		submeta['sig_sens'].append(sig_sens_stop - sig_sens_start)
	
	
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
				'power_deviation':power_deltas.real,
				'sigma':sigma,
				'start_frequency': freqs[0]*10**6,
				'middle_frequency':middlefreq*10**6,
				'axion_frequencies':axion_rmfs #in Hz
				}
	return results

        
    
    
    
    
    
    