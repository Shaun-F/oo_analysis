"""
grand_analysis.py: generates grand spectraum for a number of significant figures
Created By: Erik Lentz
"""
from oo_analysis.analysis.coaddition import add_subtract_scan, scan_cl
import h5py
from oo_analysis.analysis.scan_analysis import analysis
import numpy as np
import os

# add scan results to grand spectra
	
class analyzer(object):
	"""
	Class maintains the grand spectra. It initializes it if it doesnt exist, analysis the data, 
	adds in scans to the grand spectra, and returns the grand spectra as an output of the Grand_Analysis method
	"""
	def __init__(self, object):
		#Set class attribute that holds the HDF5 file
		self.file = object.h5py_file
		
		#Set keyword arguments
		for key, attr in object.__dict__.items():
			setattr(self, key, attr)	
		
		#pull grand spectra or create it. Then pull data group
		try:
			if "grand_spectra_run1a" in self.file.keys():
				self.grand_spectra_group = self.file["grand_spectra_run1a"]
			else:
				self.grand_spectra_group = self.file.create_group("grand_spectra_run1a")
				from oo_analysis.toolbox.grand_spectra_init import initialization
				initialization(self.grand_spectra_group, self.keys) #Initialize grand spectra datasets
		except KeyError as error:
			self.file.close()
			open(os.getcwd() + '/oo_analysis/meta/error_log', 'a+').write("\n\n"+ str(error))
			return error
			
		#Set class attribute that holds the core analysis class
		self.core_analysis = object	

		
	
	def Grand_Analysis(self):
		# cycle over scans, analyzing each
		self.analysis_results = {}
	
		#Initialize coaddition class
		grand_analysis_class = scan_cl('Null', {'Null': 'Null'}, chunk = self.grand_spectra_group, op='Null')
		
		#begin iteration
		N_iter = len(self.keys)
		### Pull keys from this class into a dictionary to be passed to the analysis
		scanparam_keys = ['restype', 'notes', 'nbins', 'axion_frequencies_to_inject',
					'pec_vel', 'signal_dataset', 'filter', 'filter_params', 'signal']

		params = dict((key, getattr(self, key)) for key in scanparam_keys)

		#run single analysis routine
		print("\rPerforming single scan analysis       ")
		self.analysis_results = dict((key, analysis(self.dig_dataset[key], scan_number=key, Tsys=self.Tsys[key],**params)) for key in self.keys if not self.dig_dataset[key].attrs['cut'])
		self.analyzed_keys = list(self.analysis_results.keys()) #List of analyzed keys
		iterator = range(len(self.analyzed_keys)) #Generator to be used in iterations

		#cycle over scans and make cuts on data

		int_times = dict((key, float(self.dig_dataset[key].attrs['integration_time'])) for key in self.analyzed_keys) #Collection of scan integration times (seconds)
		bandwidths = dict((key, 10**6*float(self.dig_dataset[key].attrs['frequency_resolution'])) for key in self.analyzed_keys) #Collection of scan resolutions (Hz)

]
		#radiometer dispersion: Tsys/root(B*t) B==bandidth t==integration time
		radiometer_dispersion = dict((key,(1/(int_times[key]*bandwidths[key])**(0.5))) for key in self.analyzed_keys) #Get collection of radiometer dispersion for each analyzed scan
		scan_dispersion = dict((key, self.analysis_results[key]['sigma']) for key in self.analyzed_keys) #Get collection of analyzed scans dispersion

		disp_toobig = dict((iter, scan_dispersion[iter]>3*radiometer_dispersion[iter]) for iter in self.analyzed_keys) #Booleans for scans whose dispersions are too large
		disp_toosmall = dict((iter, scan_dispersion[iter]<0.3*radiometer_dispersion[iter]) for iter in self.analyzed_keys) #Booleans for scans whose dispersions are too small
		cut = dict((key, bool(disp_toobig[key]|disp_toosmall[key])) for key in self.analyzed_keys) #Combined above two boolean arrays to one master array that says what dispersion are outside the bounds
		add_or_subtract= {"True": 'subtract', 
							"False":'add'} #Dictionary that holds the keywords to be passed to the coaddition scripts
		
		#Perform coaddition
		add_sub_gen = list((add_subtract_scan(add_or_subtract[str(cut[iter])], self.analysis_results[iter], grand_analysis_class, iter, self.grand_spectra_group) for iter in self.analyzed_keys))
	
		
		#Candidate flagging
		fit_significance = self.grand_spectra_group['axion_fit_significance'][...]
		axion_frequencies = self.grand_spectra_group['axion_frequencies'][...]
		candidate_mask = (fit_significance>3) 											#3-sigma cutoff
		candidates = axion_frequencies[candidate_mask] #Axion masses of candidates
		
		#Background sizes
		bg_sizes = list(self.analysis_results[key]['background_size'] for key in self.analyzed_keys)
		
		#Return the entire grand spectra group and the number of scans cut
		grand_analysis_results = (self.grand_spectra_group,
									candidates,
									len(np.where(list(cut.values()))[0]),
									bg_sizes,
									self.analyzed_keys
									)
		return grand_analysis_results
		
		
		
		
		
		