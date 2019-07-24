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
	
	def __init__(self, object):
		self.file = object.h5py_file
		#pull grand spectra or create it. Then pull data group
		try:
			if "grand_spectra_run1a" in self.file.keys():
				self.grand_spectra_group = self.file["grand_spectra_run1a"]
			else:
				self.grand_spectra_group = self.file.create_group("grand_spectra_run1a")
				from oo_analysis.toolbox.grand_spectra_init import initialization
				initialization(self.grand_spectra_group) #Initialize grand spectra datasets
		except KeyError as error:
			self.file.close()
			open(os.getcwd() + '/oo_analysis/meta/error_log', 'a+').write("\n\n"+ str(error))
			return error
		self.core_analysis = object	

		#Set keyword arguments
		for key, attr in object.__dict__.items():
			setattr(self, key, attr)	
	
	def Grand_Analysis(self):
		# cycle over scans, analyzing each
		self.analysis_results = {}
	
		grand_analysis_class = scan_cl('Null', {'Null': 'Null'}, chunk = self.grand_spectra_group, op='Null')
		
		#begin iteration
		N_iter = len(self.keys)
		scanparam_keys = ['restype', 'notes', 'nbins', 'axion_frequencies_to_inject',
					'pec_vel', 'signal_dataset', 'filter', 'filter_params', 'signal']

		params = dict((key, getattr(self, key)) for key in scanparam_keys)

		#run single analysis routine
		print("\rPerforming single scan analysis       ")
		self.analysis_results = dict((key, analysis(self.dig_dataset[key], scan_number=key, Tsys=self.Tsys[key],**params)) for key in self.keys if not self.dig_dataset[key].attrs['cut'])
		self.analyzed_keys = list(self.analysis_results.keys())
		iterator = range(len(self.analyzed_keys))

		#cycle over scans and make cuts on data

		int_times = dict((key, float(self.dig_dataset[key].attrs['integration_time'])) for key in self.analyzed_keys)
		bandwidths = dict((key, float(self.dig_dataset[key].attrs['frequency_resolution'])) for key in self.analyzed_keys)

		#radiometer dispersion: Tsys/root(B*t) B==bandidth t==integration time
		radiometer_dispersion = dict((key,(1/(int_times[key]*bandwidths[key])**(0.5))) for key in self.analyzed_keys)
		scan_dispersion = dict((key, self.analysis_results[key]['sigma']) for key in self.analyzed_keys)

		toobig_reason = {"True":"Dispersion compared to radiometer dispersion. Too large", 
						"False":""}
		toosmall_reason = {"True": "Dispersion compared to radiometer dispersion. Too small", 
							"False":""}

		disp_toobig = dict((iter, scan_dispersion[iter]>3*radiometer_dispersion[iter]) for iter in self.analyzed_keys)
		disp_toosmall = dict((iter, scan_dispersion[iter]<0.3*radiometer_dispersion[iter]) for iter in self.analyzed_keys)
		cut = dict((key, bool(disp_toobig[key]*disp_toosmall[key])) for key in self.analyzed_keys)
		add_or_subtract= {"True": 'subtract', 
							"False":'add'}
		
		add_sub_gen = (add_subtract_scan(add_or_subtract[str(cut[iter])], self.analysis_results[iter], grand_analysis_class, iter, self.grand_spectra_group) for iter in self.analyzed_keys)
		for i in add_sub_gen: #Iterate over the above generator
			pass
		
		#Return the entire grand spectra group and the number of scans cut
		return self.grand_spectra_group, len(np.where(list(cut.values()))[0])
		










