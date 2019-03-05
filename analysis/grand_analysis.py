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
import time
# add scan results to grand spectra





	
class analyser(object):
	
	def __init__(self, object):
		self.file = object.h5py_file
		#pull grand spectra or create it. Then pull data group
		try:
			grand_spectra_start = time.time()
			if "grand_spectra_run1a" in self.file.keys():
				self.grand_spectra_group = self.file["grand_spectra_run1a"]
			else:
				#self.grand_spectra = self.file.create_dataset("grand_spectra_run1a", data = [], dtype=float, chunks = True, maxshape=None)
				self.grand_spectra_group = self.file.create_group("grand_spectra_run1a")
				from toolbox.grand_spectra_init import initialization
				initialization(self.grand_spectra_group)
			grand_spectra_stop = time.time()
		except KeyError as error:
			self.file.close()
			return error
			
		if object.meta_analysis[0]:
			object.meta_analysis.append("Pulling/creating grand spectra group took {0:03f} seconds".format(grand_spectra_stop - grand_spectra_start))

		for key, attr in object.__dict__.items():
			setattr(self, key, attr)
	
		self.ncut = 0 #Counter for number of cut scans
	def Grand_Analysis(self):
		# cycle over scans, analyzing each
		self.analysis_results = {}
		
		#meta analysis
		submeta={"timeit":self.meta_analysis[0]}
		if self.meta_analysis[0]:
			meta_stats = ["constants", "axion_power", "modulation", "BS", "consolidation", "cavity_lorentz", "bin_stats", "convolutions", "axion_rmfs", "max_likelihood_arith", "max_likelihood", "sig_sens", 'reflecting_time', 'calculating_tophat_size', 'reciprocating_array', 'generating_large_array', 'generating_tophat', 'fft_and_highpassfilter', 'ifft', 'picking_og_signal', 'dividing_structure']
			submeta = {key: [] for key in meta_stats}
			submeta['timeit'] = self.meta_analysis[0]
		iteration_start = time.time()
		analysis_timer = []
		
		#begin iteration
		for key in self.keys:
			#Run single scan analysis
			scan = self.dig_dataset[key]
			scanparam_keys = ['restype', 'notes', 'nbins']
			modparam_keys = ['pec_vel', 'signal', 'filter_params']
			
			scanparams = {key: getattr(self, key) for key in scanparam_keys}
			modparams = {key: getattr(self, key) for key in modparam_keys}
			Tsys = {'Tsys': self.Tsys[key]}
			scan_num = {'scan_number':key}
			axion_scan = {'axion_scan': self.dig_dataset[key]}
			
			params = {**scanparams, **modparams, **Tsys, **scan_num, **axion_scan, "submeta":submeta}
			
			#Run single analysis routine
			analysis_start = time.time()
			times = {}
			self.analysis_results[key] =  analysis(scan, **params)
			analysis_stop = time.time()
			
			analysis_timer.append(analysis_stop - analysis_start)
			#cycle over scans and Make cuts on data
			#if deltas dispersion is factor of 3 larger than radiometer equation, cut it.
			#Radiometer dispersion: Tsys/root(B*t) B==bandwidth, t==integration time
			scan = self.dig_dataset[key]
			int_time = 100.0026931 #second             #should probably pull this from scan parameters
			bandwidth = float(scan.attrs["frequency_resolution"])*10**6 #Hz
			Tsys = self.Tsys[key]
			
			radiometer_dispersion = 1/((bandwidth*int_time)**(0.5))
			scan_dispersion = self.analysis_results[key]["sigma"]
			
			#print("\n\nRadiometer: ", radiometer_dispersion, "\nScan:", scan_dispersion, "\n fractional: ", scan_dispersion/radiometer_dispersion)
			
			cut = False
			cut_reason = ""
			if scan_dispersion>3*radiometer_dispersion:
				cut = True
				cut_reason = "Dispersion compared to radiometer dispersion. Too large"
			if scan_dispersion<0.3*radiometer_dispersion:
				cut = True
				cut_reason = "Dispersion compared to radiometer dispersion. Too Small"
				print("Background subtraction failed for scan {0}".format(key))
			scan.attrs["cut"] = cut
			scan.attrs["cut_reason"] = cut_reason
			
		
			#add remaining scan to grand_spectra via coaddition
			if scan.attrs['cut'] == False:
				add_subtract_scan('add', self.analysis_results[key], self.grand_spectra_group, key)
			elif scan.attrs['cut'] == True:
				self.ncut += 1
				add_subtract_scan('subtract', self.analysis_results[key], self.grand_spectra_group, key)
		iteration_stop = time.time()
		
		if self.meta_analysis[0]:
			avg_constants = sum(submeta['constants'])/len(submeta['constants'])
			avg_axion_power = sum(submeta['axion_power'])/len(submeta['axion_power'])
			avg_modulation = sum(submeta['modulation'])/len(submeta['modulation'])
			avg_BS = sum(submeta['BS'])/len(submeta['BS'])
			avg_consolidation = sum(submeta['consolidation'])/len(submeta['consolidation'])
			avg_cavity_lorentz = sum(submeta['cavity_lorentz'])/len(submeta['cavity_lorentz'])
			avg_bin_stats = sum(submeta['bin_stats'])/len(submeta['bin_stats'])
			avg_convolutions = sum(submeta['convolutions'])/len(submeta['convolutions'])
			avg_axion_rmfs = sum(submeta['axion_rmfs'])/len(submeta['axion_rmfs'])
			avg_max_likelihood_arith = sum(submeta['max_likelihood_arith'])/len(submeta['max_likelihood_arith'])
			avg_max_likelihood = sum(submeta['max_likelihood'])/len(submeta['max_likelihood'])
			avg_sig_sens = sum(submeta['sig_sens'])/len(submeta['sig_sens'])
			
			total_constants = sum(submeta['constants'])
			total_axion_power = sum(submeta['axion_power'])
			total_modulation = sum(submeta['modulation'])
			total_BS = sum(submeta['BS'])
			total_consolidation = sum(submeta['consolidation'])
			total_cavity_lorentz = sum(submeta['cavity_lorentz'])
			total_bin_stats = sum(submeta['bin_stats'])
			total_convolutions = sum(submeta['convolutions'])
			total_axion_rmfs = sum(submeta['axion_rmfs'])
			total_max_likelihood_arith = sum(submeta['max_likelihood_arith'])
			total_max_likelihood = sum(submeta['max_likelihood'])
			total_sig_sens = sum(submeta['sig_sens'])
			
			iteration_time = iteration_stop - iteration_start
			total_analysis_time = sum(analysis_timer)/len(analysis_timer)
			
			s0 = "\n\n########## averages of single scan analysis ##########"
			s1 = "average time to calculate constants in single scan analysis is {0:03f}".format(avg_constants)
			s2 = "average time to calculate axion power in single scan analysis is {0:03f}".format(avg_axion_power)
			s3 = "average time to perform signal modulation in single scan analysis is {0:03f}".format(avg_modulation)
			s4 = "average time to perform background subtraction in single scan analysis is {0:03f}".format(avg_BS)
			s5 = "average time to consolidate bin resolution in single scan analysis is {0:03f}".format(avg_consolidation)
			s6 = "average time to calculate cavity lorentzian in single scan analysis is {0:03f}".format(avg_cavity_lorentz)
			s7 = "average time to calculate single bin statistics in single scan analysis is {0:03f}".format(avg_bin_stats)
			s8 = "average time to calculate convolutions in single scan analysis is {0:03f}".format(avg_convolutions)
			s9 = "average time to calculate all axion frequencies in single scan analysis is {0:03f}".format(avg_axion_rmfs)
			s10 = "average time to carry out maximum likelihood arithmetic (sums, multiplications) in single scan analysis is {0:03f}".format(avg_max_likelihood_arith)
			s11 = "average time to calculate maximum likelihood in single scan analysis is {0:03f}".format(avg_max_likelihood)
			s12 = "average time to calculate sensitivities and significances in single scan analysis is {0:03f}".format(avg_sig_sens)
			
			s13 = "\n\n########## Totals of single scan analysis ##########"
			s14 = "total time to calculate constants in single scan analysis is {0:03f}".format(total_constants)
			s15 = "total time to calculate axion power in single scan analysis is {0:03f}".format(total_axion_power)
			s16 = "total time to perform signal modulation in single scan analysis is {0:03f}".format(total_modulation)
			s17 = "total time to perform background subtraction in single scan analysis is {0:03f}".format(total_BS)
			s18 = "total time to consolidate bin resolution in single scan analysis is {0:03f}".format(total_consolidation)
			s19 = "total time to calculate cavity lorentzian in single scan analysis is {0:03f}".format(total_cavity_lorentz)
			s20 = "total time to calculate single bin statistics in single scan analysis is {0:03f}".format(total_bin_stats)
			s21 = "total time to calculate convolutions in single scan analysis is {0:03f}".format(total_convolutions)
			s22 = "total time to calculate all axion frequencies in single scan analysis is {0:03f}".format(total_axion_rmfs)
			s23 = "total time to carry out maximum likelihood arithmetic (sums, multiplications) in single scan analysis is {0:03f}".format(total_max_likelihood_arith)
			s24 = "total time to calculate maximum likelihood in single scan analysis is {0:03f}".format(total_max_likelihood)
			s25 = "total time to calculate sensitivities and significances in single scan analysis is {0:03f}".format(total_sig_sens)
			
			strings = [s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25]
			for i in strings:
				self.meta_analysis.append(i)
			
			
			avg_reflecting_time = sum(submeta['reflecting_time'])/len(submeta['reflecting_time'])
			avg_calculating_tophat_size = sum(submeta['calculating_tophat_size'])/len(submeta['calculating_tophat_size'])
			avg_reciprocating_array = sum(submeta['reciprocating_array'])/len(submeta['reciprocating_array'])
			avg_generating_large_array = sum(submeta['generating_large_array'])/len(submeta['generating_large_array'])
			avg_generating_tophat = sum(submeta['generating_tophat'])/len(submeta['generating_tophat'])
			avg_fft_and_highpassfilter = sum(submeta['fft_and_highpassfilter'])/len(submeta['fft_and_highpassfilter'])
			avg_ifft = sum(submeta['ifft'])/len(submeta['ifft'])
			avg_picking_og_signal = sum(submeta['picking_og_signal'])/len(submeta['picking_og_signal'])
			avg_dividing_structure = sum(submeta['dividing_structure'])/len(submeta['dividing_structure'])
			self.meta_analysis.append("\n\n########## Averages of background subtraction ##########")
			self.meta_analysis.append("average time to reflect scan is {0:03f}".format(avg_reflecting_time))
			self.meta_analysis.append("average time to calculate tophat size is {0:03f}".format(avg_calculating_tophat_size))
			self.meta_analysis.append("average time to calculated reciprocated signal is {0:03f}".format(avg_reciprocating_array))
			self.meta_analysis.append("average time to generate large array of clones is {0:03f}".format(avg_generating_large_array))
			self.meta_analysis.append("average time to generate tophat is {0:03f}".format(avg_generating_tophat))
			self.meta_analysis.append("average time to fft and multiply by tophat is {0:03f}".format(avg_fft_and_highpassfilter))
			self.meta_analysis.append("average time to inverse fft is {0:03f}".format(avg_ifft))
			self.meta_analysis.append("average time to extract original signal is {0:03f}".format(avg_picking_og_signal))
			self.meta_analysis.append("average time to divide out structure is {0:03f}".format(avg_dividing_structure))
			self.meta_analysis.append("\n\n########## Totals of background subtraction ##########")
			total_reflecting_time = sum(submeta['reflecting_time'])
			total_calculating_tophat_size = sum(submeta['calculating_tophat_size'])
			total_reciprocating_array = sum(submeta['reciprocating_array'])
			total_generating_large_array = sum(submeta['generating_large_array'])
			total_generating_tophat = sum(submeta['generating_tophat'])
			total_fft_and_highpassfilter = sum(submeta['fft_and_highpassfilter'])
			total_ifft = sum(submeta['ifft'])
			total_picking_og_signal = sum(submeta['picking_og_signal'])
			total_dividing_structure = sum(submeta['dividing_structure'])
			self.meta_analysis.append("total time to reflect scan is {0:03f}".format(total_reflecting_time))
			self.meta_analysis.append("total time to calculate tophat size is {0:03f}".format(total_calculating_tophat_size))
			self.meta_analysis.append("total time to calculated reciprocated signal is {0:03f}".format(total_reciprocating_array))
			self.meta_analysis.append("total time to generate large array of clones is {0:03f}".format(total_generating_large_array))
			self.meta_analysis.append("total time to generate tophat is {0:03f}".format(total_generating_tophat))
			self.meta_analysis.append("total time to fft and multiply by tophat is {0:03f}".format(total_fft_and_highpassfilter))
			self.meta_analysis.append("total time to inverse fft is {0:03f}".format(total_ifft))
			self.meta_analysis.append("total time to extract original signal is {0:03f}".format(total_picking_og_signal))
			self.meta_analysis.append("total time to divide out structure is {0:03f}".format(total_dividing_structure))
			
			
		return self.grand_spectra_group, self.ncut
		










