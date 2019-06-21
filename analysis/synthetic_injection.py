#Scrip injects synthetic axion signals into data
import numpy as np
import sys 
import astropy.constans as C
import astropy.units as U
sys.path.insert(0, os.path.abspath("../")
sys.path.append("../signals")

from signal_lib import signal as sig_lib
from modulation import modulation as mod

#predefined signals to inject into data


class injector():
	
	
	def __init__(self, core_object, **kwargs):
		
		#set class definitions
		for name, arg in kwargs:
			setattr(self, name, arg)
			
		self.signal_class = sig_lib()
		self.signal_strength = 0.05/3 # Level to be recoverable with ~20 scans at SNR=3
		self.expected_dfsz_snr = 3
		self.core_object = core_object()
		
	def generate_signals(self):
		axion_freq = self.axion_frequencies
		signal_name = self.signal_name
		analysis_scan_fstarts = list(self.core_object.fstart.values())
		analysis_scan_fstop = list(self.core_object.fstop.values())
		signal_collections = [[]]*len(axion_freq)
		
		for inx, afreq in enumerate(axion_freq):
			axion_in_range_mask = np.logical_and(afreq>analysis_scan_fstarts, afreq<analysis_scan_fstop)
			keys = self.core_object.keys[axion_in_range_mask]
			datasets_toinject = {key: self.core_object.dig_dataset[key] for key in keys}
			timestamps_toinject = {key: self.core_object.timestamp[key] for key in keys}
			
			amass_ev = (afreq*U.MHz*C.hbar).to(U.ev)
			kwargs = {'pec_vel': self.pec_vel, 'timestamp':timestamps_toinject, 'signal':signal_name, 'axion_mass': amass_ev, 'dig_dataset':datasets_toinject}
			
			signal_collection[inx] = mod(**kwargs).executor()
		
		return signal_collection
















