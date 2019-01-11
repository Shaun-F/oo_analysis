"""
grand_analysis.py: generates grand spectraum for a number of significant figures

Created By: Erik Lentz
"""
import sys
sys.path.append("..")
import analysis.scan_analysis 
import analysis.coaddition
import h5py
from analysis.scan_analysis import analysis
# add scan results to grand spectra





	
class analyser(object):
	
	def __init__(self, object):
		self.file = h5py.File(b"../data/raw/run1a_data.hdf5", "r+")
		#pull grand spectra or create it. Then pull data group
		if "grand_spectra_run1a" in self.file.keys():
			self.grand_spectra = self.file["grand_spectra_run1a"]
		else:
			self.grand_spectra = self.file.create_dataset("grand_spectra_run1a", data = [], dtype=float, chunks = True, maxshape=None)

		for key, attr in object.__dict__.items():
			setattr(self, key, attr)
	
	def Grand_Analysis(self):
		# cycle over scans, analyzing each
		self.analysis_results = {}
		for key in self.keys:
			#Run single scan analysis
			scan = self.dig_dataset[key]
			scanparam_keys = ['restype', 'notes', 'nbins']
			modparam_keys = ['pec_vel', 'signal', 'filter_params']
			
			scanparams = {key: getattr(self, key) for key in scanparam_keys}
			modparams = {key: getattr(self, key) for key in modparam_keys}
			Tsys = {'Tsys': self.Tsys[key]}
			scan_num = {'scan_number':key}
			axion_scan = {'axion_scan': self.axion_dataset[key]}
			params = {**scanparams, **modparams, **Tsys, **scan_num, **axion_scan}
			
			self.analysis_results[key] =  analysis(scan, **params)
		#Make cuts on data
		#add remaining scan to grand_spectra via coaddition
		
		










