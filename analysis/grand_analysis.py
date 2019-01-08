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
	
	def __init__(self):
		self.file = h5py.File(b"../data/raw/run1a_data.hdf5", "r+")
		#pull grand spectra or create it. Then pull data group
		if "grand_spectra_run1a" in file.keys():
			self.grand_spectra = file["grand_spectra_run1a"]
		else:
			self.grand_spectra = file.create_dataset("grand_spectra_run1a", data = [], dtype=float, chunks = True, maxshape=None)

		self.data_group = object.dig_dataset
		self.data_keys = object.keys
	
	def grand_analysis(self):
		# cycle over scans, analyzing each
		self.analysis_results = {}
		for key in self.data_keys:
			#Run single scan analysis
			scan = self.data_group[key]
			scanparam_keys = ['dig_dataset', 'res', 'notes']
			modparam_keys = ['pec_vel', 'signal']
			
			scanparams = {key: getattr(object, key) for key in scanparam_keys}
			modparams = {key: getattr(object, key) for key in modparam_keys}
			params = {**scanparams, **modparams}
			
			self.analysis_results[key] =  analysis(scan, **params)
		#Make cuts on data
		
		#add remaining scan to grand_spectra via coaddition
		
		










