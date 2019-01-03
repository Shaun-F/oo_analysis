import os
import sys
from scipy.optimize import curve_fit
import numpy as np

def calc_sys_temp(h5py_digitizer_dataset, h5py_axion_dataset, squid_temp_dataset):
	"""
	Function computes the system temperature. 
	h5py_digitizer_dataset should be (h5pyfile)["digitizer_log_run1a"]["(digitizer id)"]
	h5py_axion_dataset should be (h5pyfile)["axion_log_run1a"]["(digitizer id)"]
	squid_temp_dataset should be (h5pyfile)["squid_temperature_run1a"]["timestamp"]	
	
	update 05-11-18: Internalized fetching of squid dataset
	"""
	
	dig_scan=h5py_digitizer_dataset
	axion_scan = h5py_axion_dataset
		
        #Get sensor values
	cavity_bottom_temp = float(axion_scan.attrs["cavity_top_temperature"])
	cavity_top_temp = float(axion_scan.attrs["cavity_bottom_temperature"])
	T_cavity = (cavity_bottom_temp + cavity_top_temp)/2
	timestamp = axion_scan.attrs["timestamp"]
	
	squid_temperature = float(squid_temp_dataset[...]) #This is the temperature of A4 in the RF chain.
        
	fitted_func = power_fitter(dig_scan, axion_scan)
	R = np.max(fitted_func)/np.min(fitted_func)

	System_temp = (squid_temperature - T_cavity)/(R-1) #See supplementary material of the run1a results paper for equation explanation

	return {"system temperature":System_temp, "power ratio":R, "squid temperature":squid_temperature, "cavity temperature":T_cavity}


def lorentzian(resonant_frequency, Q, size_of_lorentzian_in_bins, bin_size):
	length = size_of_lorentzian_in_bins
	bin = bin_size
	res_freq = resonant_frequency
	freqs = np.arange(start = res_freq - (length/2), stop = res_freq + (length/2), step = bin, dtype=np.dtype(float))

	lorentzian = (1/(2*np.pi))*(res_freq/Q)/((freqs-res_freq)**2 + (res_freq/(2*Q))**2)
	return lorentzian


def power_fitter(digitizer_scan, axion_scan):

	dig_scan = digitizer_scan
	power_spec = dig_scan[...] #power spectrum ch 1
	bin_size = float(dig_scan.attrs["frequency_resolution"])   #frequency res ch1. in units of Hz
	Q = float(axion_scan.attrs["Q"]) #Q ch1
	res_freq = float(axion_scan.attrs["mode_frequency"])     #mode freq ch1. in units of Hz
	length= len(power_spec)

	lorentzian = lambda f: (1/(2*np.pi))*(res_freq/Q)/((f-res_freq)**2 + (res_freq/(2*Q))**2)
	fitting_func = lambda f,a_0,a_1,a_2,a_3: (a_0 + a_1*lorentzian(f) + a_2*(f-res_freq)*lorentzian(f))*(1+a_3*(f-res_freq))

	xdata = np.arange(start = res_freq - (length/2)*bin_size, stop = res_freq + (length/2)*bin_size-bin_size, step = bin_size, dtype=np.dtype(float))
	minimized_params, covariance = curve_fit(fitting_func, xdata, power_spec, p0=[np.max(power_spec),0,0,0], )
	return fitting_func(xdata, minimized_params[0], minimized_params[1], minimized_params[2], minimized_params[3])
	
def get_squid_dataset(timestamp):
	"""function gets the squid_temperature value closest to the input timestamp. Must be in format dd-mm-yyyy hh:mm:ss"""
	f = h5py.File('run1a_data', 'r')["squid_temperature_run1a"]
	keys = [x for x in f.keys()]
	year, month, day = str(timestamp[6:10]),str(timestamp[3:5]),str(timestamp[0:2])
	hour, minute, second = str(timestamp[11:13]),str(timestamp[14:16]),str(timestamp[17:19])
	for x in keys:
		if (x[6:10]==year and x[3:5]==month and x[0:2]==day):
			if x[11:13]==hour:
				if x[14:16]==minute:
					return f[x]
				else:
					for i in np.arange(15):
						minute = int(minute) #turn into number
						if (x[14:16]==str(minute-i) or x[14:16]==str(minute+i)):
							#Find closest value
							return f[x]
