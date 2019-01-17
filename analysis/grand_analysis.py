"""
grand_analysis.py: generates grand spectraum for a number of significant figures

Created By: Erik Lentz
"""
import sys
sys.path.append("..")
import analysis.scan_analysis 
import analysis.coaddition; from analysis.coaddition import add_subtract_scan
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

			#cycle over scans and Make cuts on data
			#if deltas dispersion is factor of 3 larger than radiometer equation, cut it.
			#Radiometer dispersion: Tsys/root(B*t) B==bandwidth, t==integration time
			scan = self.dig_dataset[key]
			int_time = 100.0026931 #second             #should probably pull this from scan parameters
			bandwidth = float(scan.attrs["frequency_resolution"])*10**6 #Hz
			Tsys = self.Tsys[key]
			
			radiometer_dispersion = Tsys/((bandwidth*int_time)**(0.5))
			scan_dispersion = self.analysis_results[key]["sigma"]
			
			cut = False
			cut_reason = ""
			if scan_dispersion>3*radiometer_dispersion:
				cut = True
				cut_reason = "Dispersion compared to radiometer dispersion. To large"
				scan.attrs["cut"] = cut
				scan.attrs["cut_reason"] = cut_reason
		
			#add remaining scan to grand_spectra via coaddition
			if scan.attrs['cut'] == False:
				add_subtract_scan('add', self.analysis.results[key], self.grand_spectra, key)
			elif scan.attrs['cut'] == True:
				add_subtract_scan('subtract', self.analysis.results[key], self.grand_spectra, key)
				
		return self.grand_spectra
		










