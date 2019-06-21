"""
grand_analysis.py: generates grand spectraum for a number of significant figures

Created By: Erik Lentz
"""
import sys
sys.path.append("..")
import analysis.scan_analysis 
import analysis.coaddition; from analysis.coaddition import add_subtract_scan, scan_cl
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
			open('../meta/error_log', 'a+').write(str(time.time())+ "\n\n"+ str(error))
			return error
		self.core_analysis = object	
		if object.meta_analysis[0]:
			object.meta_analysis.append("\nPulling/creating grand spectra group took {0:03f} seconds".format(grand_spectra_stop - grand_spectra_start))

		for key, attr in object.__dict__.items():
			setattr(self, key, attr)
	
		self.ncut = 0 #Counter for number of cut scans
	
	
	def Grand_Analysis(self):
		# cycle over scans, analyzing each
		self.analysis_results = {}
		
		#meta analysis
		submeta={"timeit":self.meta_analysis[0]}
		if self.meta_analysis[0]:
			meta_stats = ["constants", "axion_power", "modulation", "BS", 'convolutions', "consolidation", "cavity_lorentz", "bin_stats", "axion_rmfs", "sig_sens", 'reflecting_time', 
							'calculating_tophat_size', 'reciprocating_array', 'generating_large_array', 'generating_tophat', 'fft_and_highpassfilter', 'ifft', 
							'picking_og_signal', 'dividing_structure', 'consolidation_time', 'class_init_time', 'optimal_weight_sum_consolidation_time', 
							'model_excess_sqrd_consolidation_time', 'nscan_consolidation_time', 'SNR_consolidation_time', 'power_deviation_consolidation_time', 
							'sigma_A_consolidation_time', 'coupling_sensitivity_consolidation_time', 'noise_power_consolidation_time', 'axion_fit_consolidation_time', 
							'axion_fit_significance_consolidation_time', 'power_sensitivity_consolidation_time',
							'init_scans_in_class', 'init_grand_spectra', 'init_grand_spectra_in_class', 'growing_grand_spectra', 'reinit_grand_spectra_in_class', 'last_calc_time',
							'attaching_new_scan_to_class','calculating_indices','close_out']
			submeta = {key: [] for key in meta_stats}
			submeta['timeit'] = self.meta_analysis[0]
		
		iteration_start = time.time()
		analysis_timer = []
		coaddition_timer = []
		
		grand_analysis_class = scan_cl('Null', {'Null': 'Null'}, chunk = self.grand_spectra_group, op='Null', **{'submeta':submeta})
		
		#begin iteration
		for key in self.keys:
			try:
				#Run single scan analysis
				scan = self.dig_dataset[key]
				scanparam_keys = ['restype', 'notes', 'nbins']
				modparam_keys = ['pec_vel', 'signal_dataset', 'filter_params']
				
				scanparams = {key: getattr(self, key) for key in scanparam_keys}
				modparams = {key: getattr(self, key) for key in modparam_keys}
				Tsys = {'Tsys': self.Tsys[key]}
				scan_num = {'scan_number':key}
				axion_scan = {'axion_scan': self.dig_dataset[key]}
				
				params = {**scanparams, **modparams, **Tsys, **scan_num, **axion_scan, "submeta":submeta}
				
				#Run single analysis routine
				analysis_start = time.time()
				self.analysis_results[key] =  analysis(scan, **params)
				analysis_stop = time.time()
				
				analysis_timer.append(analysis_stop - analysis_start)
				#cycle over scans and Make cuts on data
				#if deltas dispersion is factor of 3 larger than radiometer equation, cut it.
				#Radiometer dispersion: Tsys/root(B*t) B==bandwidth, t==integration time
				scan = self.dig_dataset[key]
				int_time = float(scan.attrs['integration_time']) #second
				bandwidth = float(scan.attrs["frequency_resolution"])*10**6 #Hz
				
				radiometer_dispersion = 1/((bandwidth*int_time)**(0.5))
				scan_dispersion = self.analysis_results[key]["sigma"]
				
				#print("\n\nRadiometer: ", radiometer_dispersion, "\nScan:", scan_dispersion, "\n fractional: ", scan_dispersion/radiometer_dispersion)
				cut = False
				cut_reason = ""
				if scan_dispersion>3*radiometer_dispersion:
					cut = True
					cut_reason = "Dispersion compared to radiometer dispersion. Too large"
					self.core_analysis.bad_scans.append(key)
				if scan_dispersion<0.3*radiometer_dispersion:
					cut = True
					cut_reason = "Dispersion compared to radiometer dispersion. Too Small"
					self.core_analysis.bad_scans.append(key)
					print("Background subtraction failed for scan {0}".format(key))
				
				if not scan.attrs['cut']:
					scan.attrs["cut"] = cut
					scan.attrs["cut_reason"] = cut_reason
					
				
				#add remaining scans to grand_spectra via coaddition
				coaddition_start = time.time()
				if not scan.attrs['cut']:
					add_subtract_scan('add', self.analysis_results[key], grand_analysis_class, key, self.grand_spectra_group, **{'submeta':submeta})
				
				elif scan.attrs['cut'] and scan.attrs['cut_reason']!="scan outside 645 to 685 MHz range":
					self.ncut += 1
					add_subtract_scan('subtract', self.analysis_results[key], grand_analysis_class, key, self.grand_spectra_group, **{'submeta':submeta})
				coaddition_stop = time.time()
				coaddition_timer.append(coaddition_stop - coaddition_start)
				
			except (KeyError, MemoryError, IndexError) as error:
				print("Error at scan {0}. Saving to error log".format(key))
				open('../meta/error_log', 'a+').write(str(time.time())+ "\n\n"+ str(error))
				self.file.close()
				raise
		iteration_stop = time.time()
		
		if self.meta_analysis[0]:
			def avg(list):
				if len(list)==0:
					return 0
				return sum(list)/len(list)
			avg_constants = sum(submeta['constants'])/len(submeta['constants'])
			avg_axion_power = sum(submeta['axion_power'])/len(submeta['axion_power'])
			avg_modulation = sum(submeta['modulation'])/len(submeta['modulation'])
			avg_BS = sum(submeta['BS'])/len(submeta['BS'])
			avg_consolidation = sum(submeta['consolidation'])/len(submeta['consolidation'])
			avg_cavity_lorentz = sum(submeta['cavity_lorentz'])/len(submeta['cavity_lorentz'])
			avg_bin_stats = sum(submeta['bin_stats'])/len(submeta['bin_stats'])
			avg_convolutions = sum(submeta['convolutions'])/len(submeta['convolutions'])
			avg_axion_rmfs = sum(submeta['axion_rmfs'])/len(submeta['axion_rmfs'])
			
						
			avg_coaddition = avg(coaddition_timer)
			avg_consolidation_time = avg(submeta['consolidation_time'])
			#avg_scans_in_time = avg(submeta['scans_in_time'])
			#avg_class_init_time = avg(submeta['class_init_time'])
			avg_opt_wght_sum_cons_time = avg(submeta['optimal_weight_sum_consolidation_time'])
			avg_model_exc_sqrd_cons_time = avg(submeta['model_excess_sqrd_consolidation_time'])
			avg_nscan_cons_time = avg(submeta['nscan_consolidation_time'])
			avg_SNR_cons_time = avg(submeta['SNR_consolidation_time'])
			avg_sig_A_cons_time = avg(submeta['sigma_A_consolidation_time'])
			avg_pow_dev_cons_time = avg(submeta['power_deviation_consolidation_time'])
			avg_coup_sens_cons_time = avg(submeta['coupling_sensitivity_consolidation_time'])
			avg_pow_sens_cons_time = avg(submeta['power_sensitivity_consolidation_time'])
			avg_noise_pow_cons_time = avg(submeta['noise_power_consolidation_time'])
			avg_axion_fit_cons_time = avg(submeta['axion_fit_consolidation_time'])
			avg_axion_fit_sig_cons_time = avg(submeta['axion_fit_significance_consolidation_time'])
			#avg_scans_in_GS_add_time = avg(submeta['scans_in_grand_spectra_addition_time'])
			avg_last_calc_time = avg(submeta['last_calc_time'])
			avg_attaching_new_scan_time = avg(submeta['attaching_new_scan_to_class'])
			avg_calculating_indices_time = avg(submeta['calculating_indices'])
			#avg_check_member_time = avg(submeta['check_member_time'])
			
			avg_init_scans_in_class_time = avg(submeta['init_scans_in_class'])
			avg_init_GS = avg(submeta['init_grand_spectra'])
			avg_init_GS_in_cl = avg(submeta['init_grand_spectra_in_class'])
			avg_close_out_time = avg(submeta['close_out'])
			#avg_growing_GS = avg(submeta['growing_grand_spectra'])
			#avg_reinit_GS_in_cl = avg(submeta['reinit_grand_spectra_in_class'])
			
			total_constants = sum(submeta['constants'])
			total_axion_power = sum(submeta['axion_power'])
			total_modulation = sum(submeta['modulation'])
			total_BS = sum(submeta['BS'])
			total_consolidation = sum(submeta['consolidation'])
			total_cavity_lorentz = sum(submeta['cavity_lorentz'])
			total_bin_stats = sum(submeta['bin_stats'])
			total_convolutions = sum(submeta['convolutions'])
			total_axion_rmfs = sum(submeta['axion_rmfs'])
			
			total_coaddition = sum(coaddition_timer)
			total_consolidation_time = sum(submeta['consolidation_time'])
			#total_scans_in_time = sum(submeta['scans_in_time'])
			#total_class_init_time = sum(submeta['class_init_time'])
			total_opt_wght_sum_cons_time = sum(submeta['optimal_weight_sum_consolidation_time'])
			total_model_exc_sqrd_cons_time = sum(submeta['model_excess_sqrd_consolidation_time'])
			total_nscan_cons_time = sum(submeta['nscan_consolidation_time'])
			total_SNR_cons_time = sum(submeta['SNR_consolidation_time'])
			total_sig_A_cons_time = sum(submeta['sigma_A_consolidation_time'])
			total_pow_dev_cons_time = sum(submeta['power_deviation_consolidation_time'])
			total_coup_sens_cons_time = sum(submeta['coupling_sensitivity_consolidation_time'])
			total_pow_sens_cons_time = sum(submeta['power_sensitivity_consolidation_time'])
			total_noise_pow_cons_time = sum(submeta['noise_power_consolidation_time'])
			total_axion_fit_cons_time = sum(submeta['axion_fit_consolidation_time'])
			total_axion_fit_sig_cons_time = sum(submeta['axion_fit_significance_consolidation_time'])
			#total_scans_in_GS_add_time = sum(submeta['scans_in_grand_spectra_addition_time'])
			total_last_calc_time = sum(submeta['last_calc_time'])
			total_attaching_new_scan_time = sum(submeta['attaching_new_scan_to_class'])
			total_calculating_indices_time = sum(submeta['calculating_indices'])
			total_close_out_time = sum(submeta['close_out'])
			#total_check_member_time = sum(submeta['check_member_time'])
			
			total_init_scans_in_class_time = sum(submeta['init_scans_in_class'])
			total_init_GS = sum(submeta['init_grand_spectra'])
			total_init_GS_in_cl = sum(submeta['init_grand_spectra_in_class'])
			total_growing_GS = sum(submeta['growing_grand_spectra'])
			total_reinit_GS_in_cl = sum(submeta['reinit_grand_spectra_in_class'])
			
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
			
			
			s13 = "average time to add/subtract scan is {0:03f}".format(avg_coaddition)
			s14 = "average time to consolidation all data is {0:03f}".format(avg_consolidation_time)
			s15 = ''#"average time to add/subtract scan id is {0:03f}".format(avg_scans_in_time)
			s16 = ''#"average time to initialize coaddition class is {0:03f}".format(avg_class_init_time)
			s17 = "average time to consolidate optimal_weight_sum statistic is {0:03f}".format(avg_opt_wght_sum_cons_time)
			s18 = "average time to consolidate model_excess_squared statistic is {0:03f}".format(avg_model_exc_sqrd_cons_time)
			s19 = "average time to consolidate nscan statistic is {0:03f}".format(avg_nscan_cons_time)
			s20 = "average time to consolidate SNR statistic is {0:03f}".format(avg_SNR_cons_time)
			s21 = "average time to consolidate sigma_A statistic is {0:03f}".format(avg_sig_A_cons_time)
			s22 = "average time to consolidate power_deviation statistic is {0:03f}".format(avg_pow_dev_cons_time)
			s23 = "average time to consolidate coupling_sensitivity statistic is {0:03f}".format(avg_coup_sens_cons_time)
			s24 = "average time to consolidate power_sensitivity statistic is {0:03f}".format(avg_pow_sens_cons_time)
			s25 = "average time to consolidate noise_power statistic is {0:03f}".format(avg_noise_pow_cons_time)
			s26 = "average time to consolidate axion_fit statistic is {0:03f}".format(avg_axion_fit_cons_time)
			s27 = "average time to consolidate axion_fit_significance statistic is {0:03f}".format(avg_axion_fit_sig_cons_time)
			s28 = ''#"average time to add/subtract scans_in_grand_spectra id is {0:03f}".format(avg_scans_in_GS_add_time)
			a = "average time to initialize scans in coaddition class is {0:03f}".format(avg_init_scans_in_class_time)
			b = "average time to initialize grand_spectra is {0:03f}".format(avg_init_GS)
			c = "average time to initialize grand spectra in coaddition class is {0:03f}".format(avg_init_GS_in_cl)
			d = "" #"average time to grow grand spectra is {0:03f}".format(avg_growing_GS
			e = "average time to input last_calculation time is {0:0.3f}".format(avg_last_calc_time)
			f = 'average time to calculate position of scan in grand spectra is {0:0.3f}'.format(avg_calculating_indices_time)
			g = 'average time to attach new scan to coaddition class is {0:0.3f}'.format(avg_attaching_new_scan_time)
			h = 'average time to save grand spectra file is {0:0.3f}'.format(avg_close_out_time)
			
			
			s29 = "\n\n########## Totals of single scan analysis ##########"
			s30 = "total time to calculate constants in single scan analysis is {0:03f}".format(total_constants)
			s31 = "total time to calculate axion power in single scan analysis is {0:03f}".format(total_axion_power)
			s32	= "total time to perform signal modulation in single scan analysis is {0:03f}".format(total_modulation)
			s33 = "total time to perform background subtraction in single scan analysis is {0:03f}".format(total_BS)
			s34 = "total time to consolidate bin resolution in single scan analysis is {0:03f}".format(total_consolidation)
			s35 = "total time to calculate cavity lorentzian in single scan analysis is {0:03f}".format(total_cavity_lorentz)
			s36 = "total time to calculate single bin statistics in single scan analysis is {0:03f}".format(total_bin_stats)
			s37 = "total time to calculate convolutions in single scan analysis is {0:03f}".format(total_convolutions)
			s38 = "total time to calculate all axion frequencies in single scan analysis is {0:03f}".format(total_axion_rmfs)
			
			
			
			s42 = "total time to add/subtract scan is {0:03f}".format(total_coaddition)
			s43	= "total time to consolidation all data is {0:03f}".format(total_consolidation_time)
			s44 = ''#"total time to add/subtract scan id is {0:03f}".format(total_scans_in_time)
			s45 = ''#"total time to initialize coaddition class is {0:03f}".format(total_class_init_time)
			s46 = "total time to consolidate optimal_weight_sum statistic is {0:03f}".format(total_opt_wght_sum_cons_time)
			s47 = "total time to consolidate model_excess_squared statistic is {0:03f}".format(total_model_exc_sqrd_cons_time)
			s48 = "total time to consolidate nscan statistic is {0:03f}".format(total_nscan_cons_time)
			s49 = "total time to consolidate SNR statistic is {0:03f}".format(total_SNR_cons_time)
			s50 = "total time to consolidate sigma_A statistic is {0:03f}".format(total_sig_A_cons_time)
			s51 = "total time to consolidate power_deviation statistic is {0:03f}".format(total_pow_dev_cons_time)
			s52 = "total time to consolidate coupling_sensitivity statistic is {0:03f}".format(total_coup_sens_cons_time)
			s53	= "total time to consolidate power_sensitivity statistic is {0:03f}".format(total_pow_sens_cons_time)
			s54 = "total time to consolidate noise_power statistic is {0:03f}".format(total_noise_pow_cons_time)
			s55 = "total time to consolidate axion_fit statistic is {0:03f}".format(total_axion_fit_cons_time)
			s56 = "total time to consolidate axion_fit_significance statistic is {0:03f}".format(total_axion_fit_sig_cons_time)
			s57 = ''#"total time to add/subtract scans_in_grand_spectra id is {0:03f}".format(total_scans_in_GS_add_time)
			
			i = "total time to initialize scans in coaddition class is {0:03f}".format(total_init_scans_in_class_time)
			j = "total time to initialize grand_spectra is {0:03f}".format(total_init_GS)
			k = "total time to initialize grand spectra in coaddition class is {0:03f}".format(total_init_GS_in_cl)
			l = "total time to input last_calculation time is {0:0.3f}".format(total_last_calc_time)
			m = 'total time to calculate position of scan in grand spectra is {0:0.3f}'.format(total_calculating_indices_time)
			n = 'total time to attach new scan to coaddition class is {0:0.3f}'.format(total_attaching_new_scan_time)
			o = 'total time to save grand spectra to file is {0:0.3f}'.format(total_close_out_time)
			
			
			strings = [s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,a,b,c,d,e,s29,s30,s31,s32,s33,s34,s35,s36,s37,s38,s42,s43,s44,s45,s46,s47,s48
			,s49,s50,s51,s52,s53,s54,s55,s56,s57, a,b,c,d,e,f,g,h,i,j,k,l,m,n,o]
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
		










