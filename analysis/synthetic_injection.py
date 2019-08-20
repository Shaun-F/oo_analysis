#Scrip injects synthetic axion signals into data
import numpy as np
import pandas as pd
import sys
import os
import astropy.constants as C
import astropy.units as U
import copy 
from oo_analysis.signals.signal_lib import signal as sig_lib
from oo_analysis.signals.modulation import modulation as mod
from oo_analysis.toolbox.lorentzian import lorentzian
from oo_analysis.toolbox.axion_power import axion_power

#predefined signals to inject into data


#### Procedure ####
# First pull in some file that holds the axion frequencies to inject
# Find the set of scans whose frequency range overlaps with the axion frequencies
# Use the timestamps of those scans, generate the signals modulated by the diurnal variation of the experimental velocity


class axion_injector():
	
	
	def __init__(self, dataset, filter_shape, **kwargs):
		
		#set class definitions
		for name, arg in kwargs.items():
			setattr(self, name, arg)
			
		self.signal_strength = 2.7 #Fraction of power relative to DFSZ
		self.dataset = dataset
		self.filter_shape = filter_shape
		
		#parse dataset and define axion frequencies
		self.axion_frequencies_hz = self.axion_frequencies_to_inject #List of axion frequencies to inject in units of hertz
		self.analysis_scan_fstart = float(self.dataset.attrs['start_frequency'])*10**6 #Get the list of starting frequencies
		self.analysis_scan_fstop = float(self.dataset.attrs['stop_frequency'])*10**6 #Get the list of ending frequencies
		self.analysis_scan_fres = float(self.dataset.attrs['frequency_resolution'])*10**6 #Get the list of frequency resolutions
		self.mode_frequency = float(self.dataset.attrs['mode_frequency'])*10**6 #Central frequency
		self.timestamp = self.dataset.attrs['timestamp'] #Timestamp of scan for signal generation
		self.scan_note = self.dataset.attrs['notes'] #Notes of scan 
		self.scan_data = copy.deepcopy(self.dataset[...])
		self.scan_data[0] = np.mean([self.scan_data[1], self.scan_data[2], self.scan_data[3]]) #Rescale beginning of scan to remove sampling errors
		self.Q = float(self.dataset.attrs['Q'])
		self.Tsys = float(self.dataset.attrs['squid_temperature'])
		self.generate_signal()
		
	def generate_signal(self):
		
		
		signal_collection = {} #Collection of signal distributions categorized by the axion rest mass frequency
		self.axion_in_range_mask = {} #Mask the datasets that contain the injected signal
		self.datasets_toinject = {} #Collection of datasets to inject an axion into
		self.timestamps_toinject = {} #Collection of the timestamps from the datasets which contain an injected axion
		self.bools = {} #Collection of booleans indicating if axion with frequency given by the dict key should be injected.
		
		for inx, afreq in enumerate(self.axion_frequencies_hz): #Iterate over the axion rest mass frequencies
			self.bools[str(afreq)] = (np.logical_and(afreq>self.analysis_scan_fstart, afreq<self.analysis_scan_fstop)) #Generate the mask to pick out the datasets that contain the axion rest mass frequency
			
			if self.bools[str(afreq)]:
				amass_ev = (afreq*U.Hz*C.hbar).to(U.eV).value
				kwargs = {'pec_vel': self.pec_vel, 
							'timestamp':{self.dataset.name[-6:]: self.timestamp}, 
							'signal':self.signal, 
							'axion_mass': amass_ev, 
							'dig_dataset':self.dataset, 
							'keys':[self.dataset.name[-6:]],
							'SIA':True,
							"meta_analysis":[False]}
				signal_collection[str(afreq)] = [i['signal'] for i in list(mod(**kwargs).executor().values())] #Generate the signal for each timestamp.
	
		self.signal_collection = signal_collection #store signals in an class instance attribute
		return None

	def inject_signal(self):
	
		for afreq, signal_set in self.signal_collection.items(): #Iterate over the collection of axion Rest mass frequencies.
			
			for signal_inx, signal in enumerate(signal_set): #Iterate over the collection of signals at this particular axion rest mass frequency
								
				if self.bools[str(afreq)] & (pd.isnull(self.scan_note) or not ('Synthetic axion injected' in self.notes)): #If axion already injected, dont do it again.
					freq_list = np.arange(self.analysis_scan_fstart, self.analysis_scan_fstop, self.analysis_scan_fres) #Frequency list for scan
					
					axion_starting_index_in_scan = np.argmin(np.abs(float(afreq)-freq_list)) #Get the index where the axion signal will begin in the scan data
					power_axion = axion_power(self.dataset, float(afreq))
					expected_dfsz_snr = (power_axion/(C.k_B*self.analysis_scan_fres*self.Tsys)).value #Expected signal to noise ratio
					
					if len(self.scan_data)-len(signal)-axion_starting_index_in_scan<0:
						cutoff = len(signal)+len(self.scan_data)-len(signal)-axion_starting_index_in_scan
						signal = np.asarray(np.append(np.zeros(axion_starting_index_in_scan), signal[:cutoff]))
					else:
						signal = np.asarray(np.append(np.append(np.zeros(axion_starting_index_in_scan), signal), np.zeros(len(self.scan_data)-len(signal)-axion_starting_index_in_scan))) #Append zeros to the signal to match up with the scan data
					
					power_to_add = self.filter_shape*self.signal_strength*expected_dfsz_snr*signal*lorentzian(self.Q, self.mode_frequency, self.analysis_scan_fstart, self.analysis_scan_fstop, self.analysis_scan_fres) #Calculate the power to inject into the dataset
					##Need to multiply by background shape
					
					self.scan_data += power_to_add #inject the power
					
					
					
					if pd.isnull(self.scan_note):
						self.dataset.attrs['notes'] = " Synthetic axion injected with frequency {0}".format(afreq)
					else:
					
						self.dataset.attrs['notes'] += " Synthetic axion injected with frequency {0}".format(afreq)
					
	
				else:
					pass
		SIA_meta = {'signal_collection':self.signal_collection, 'is injected':self.bools}
		return self.scan_data, SIA_meta
