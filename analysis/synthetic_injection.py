#Scrip injects synthetic axion signals into data
import numpy as np
import sys; import os
import astropy.constants as C
import astropy.units as U
sys.path.insert(0, os.path.abspath("../"))
sys.path.append("../signals")

from signal_lib import signal as sig_lib
from modulation import modulation as mod
from toolbox.lorentzian import lorentzian

#predefined signals to inject into data


#### Procedure ####
# First pull in some file that holds the axion frequencies to inject
# Find the set of scans whose frequency range overlaps with the axion frequencies
# Use the timestamps of those scans, generate the signals modulated by the diurnal variation of the experimental velocity


class axion_injector():
	
	
	def __init__(self, core_object, **kwargs):
		
		#set class definitions
		for name, arg in kwargs:
			setattr(self, name, arg)
			
		self.signal_strength = 0.05/3 # Level to be recoverable with ~20 scans at SNR=3
		self.expected_dfsz_snr = 3
		self.core_object = core_object
		
		axion_frequencies_file = open(os.path.dirname(os.path.realpath(__file__))+"/software_injection_frequencies.dat", "r")
		self.axion_frequencies_hz = list(map(int,axion_frequencies_file.readlines())) #List of axion frequencies to inject in units of hertz
		self.generate_signals()
		
	def generate_signals(self):
		signal_name = self.core_object.signal
		analysis_scan_fstarts = np.asarray(list(self.core_object.fstart.values()))
		analysis_scan_fstop = np.asarray(list(self.core_object.fstop.values()))
		
		signal_collection = {}
		self.axion_in_range_mask = {}
		self.datasets_toinject = {}
		self.timestamps_toinject = {}
		
		for inx, afreq in enumerate(self.axion_frequencies_hz):
			self.axion_in_range_mask[str(afreq)] = np.logical_and(afreq>analysis_scan_fstarts, afreq<analysis_scan_fstop)
			keys = np.asarray(self.core_object.keys)[self.axion_in_range_mask[str(afreq)]]
			if len(keys)!=0:
				self.datasets_toinject[str(afreq)] = np.asarray(list(self.core_object.dig_dataset.values()))[self.axion_in_range_mask[str(afreq)]]
				for key in keys:
					self.core_object.dig_dataset[key].attrs['notes']= "Software injection with axion frequency = {0:0.1}".format(float(afreq))
					
				self.timestamps_toinject[str(afreq)] = np.asarray(list(self.core_object.timestamp.values()))[self.axion_in_range_mask[str(afreq)]]
				
				print(self.timestamps_toinject)
				amass_ev = (afreq*U.MHz*C.hbar).to(U.eV).value
				kwargs = {'pec_vel': self.core_object.pec_vel, 'timestamp':self.timestamps_toinject[str(afreq)], 'signal':signal_name, 'axion_mass': amass_ev, 'dig_dataset':self.datasets_toinject[str(afreq)]}
				
				signal_collection[str(afreq)] = mod(**kwargs).executor()
		
		self.signal_collection = signal_collection #store signals in an class instance attribute
		return None

	def inject_signal(self):
		
		for afreq, signal_set in self.signal_collection.items():
			for signal_inx, signal in enumerate(signal_set):
				scan_data = np.asarray(self.datasets_toinject[str(afreq)][signal_inx][...])
				Q = self.datasets_toinject[str(afreq)][signal_inx].attrs['Q']
				mode_freq = self.datasets_toinject[str(afreq)][signal_inx].attrs['mode_frequency']
				fstart_scan = analysis_scan_fstarts[self.axion_in_range_mask[str(afreq)]][signal_inx]
				fstop_scan = analysis_scan_fstop[self.axion_in_range_mask[str(afreq)]][signal_inx]
				
				axion_index_in_scan = np.argmin(np.abs(float(afreq)-scan_data))
				signal = np.asarray(np.append(np.append(np.zeros(axion_index_in_scan), signal), np.zeros(len(scan_data-len(signal)-axion_index_in_scan))))
				
				power_to_add = self.signal_strength*self.expected_dfsz_snr*signal*lorentzian(Q, mode_freq, fstart_scan, fstop_scan, 95.4)
				scan_data +=power_to_add
				self.datasets_toinject[str(afreq)][signal_inx][...] = scan_data
				
		






