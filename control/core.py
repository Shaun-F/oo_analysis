"""
core.py: holds  control routines for analysis

Created by: Erik Lentz and Modified by Shaun Fell
Creation Date: 6/1/18
"""
import sys; import os
from oo_analysis.control.param_parser import parser
from oo_analysis.experiment.get_squid_dataset import get_squid_dataset
from oo_analysis.experiment.calc_sys_temp_offline import calc_sys_temp
import time; import datetime
import argparse
import copy 
from oo_analysis.toolbox.plot_dataset import plotter
from oo_analysis.toolbox.freq_to_mass import freq_to_mass
from oo_analysis.analysis.synthetic_injection import axion_injector
import urllib
from pandas import isnull
import numpy as np
import pandas as pd
import datetime as dt

# create main structure
class core_analysis():
	
	
	def __init__(self, args=None, **kwargs):
        # get parameters
		
		current_dir = os.getcwd()
		run_definitions_filename = current_dir + "/oo_analysis/job.param"  # from command line, nominally
		axion_frequencies_filename = current_dir + "/oo_analysis/analysis/software_injection_frequencies.dat"
		self.bad_scans_filename = current_dir + "/oo_analysis/analysis/bad_scans.dat"
		self.output_file = current_dir + "/oo_analysis/output/grand_spectra.dat"
		self.error_file = current_dir + "/oo_analysis/meta/error_log"
		
		self.params = parser(run_definitions_filename)
        # set switches of analysis
        # default
		self.make_plots = False
        # from param file
		
		#parse command line arguments
		if args!=None:
			args_dict = vars(args)
			[setattr(self, key, args_dict[key]) for key in list(args_dict.keys())]
			self.params['start_scan'] = self.start_scan
			self.params['end_scan'] = self.end_scan
		
		#add class attributes
		for arg, val in kwargs.items():
			if arg in self.params.keys():
				self.params[str(arg)] = val
		
		run_definitions_string = " Run Definitions: \n \t filter = {0} \n \t signal = {1} \n \t annual modulation = {2} \n \t resolution type = {3} \n \t starting scan = {4} \n \t ending scan = {5} \n ".format(self.params['filter'], self.params['signal'], self.params['pec_vel'], self.params['restype'], self.params['start_scan'], self.params['end_scan'])
		print(run_definitions_string)		
		
		for key,value in self.params.items(): # for python 3.x
			setattr(self,key,value)
		#find bad scans saved to file
		self.bad_scans_file = open(self.bad_scans_filename, 'r')
		self.bad_scans = list(map(int, self.bad_scans_file.readlines()))
		self.params['bad_scans'] = self.bad_scans
		
        # set data structures
		from oo_analysis.data import input, add_input
		pulldata_start = time.time()
		self.dig_dataset, self.h5py_file, self.no_axion_log, self.partitioned = input(self.params)
		self.keys = [copy.deepcopy(i) for i in self.dig_dataset.keys() if i not in self.no_axion_log] #Was originally dig_dataset,but some digitizer logs didnt have associated axion logs.
		pulldata_stop = time.time()
		
		print("Loading data successful. It took {0:0.3f} seconds. Beginning analysis of {1} scans".format((pulldata_stop-pulldata_start), len(self.keys)))
		if self.partitioned:
			print("Note: Dataset was partitioned")
		
		self.update_astronomical_tables()
        
		
		# derive necessary experiment data structures (put into dig_dataset)
		#Populate parameter dict's with dataset attributes.
		self.Tsys = {key: copy.deepcopy(self.dig_dataset[key].attrs["squid_temperature"]) for key in self.keys} #temperature of system during scan
		self.timestamp = {key: copy.deepcopy(self.dig_dataset[key].attrs["alog_timestamp"]) for key in self.keys} #timestamp of scan
		self.mode_frequencies = {key: copy.deepcopy(self.dig_dataset[key].attrs['mode_frequency']) for key in self.keys} #mode frequency of scan
		self.axion_mass = {key: freq_to_mass(copy.deepcopy(self.mode_frequencies[key])*10**6) for key in self.keys} #axion mass in eV with corresponding frequency equal to mode frequency
		self.fstart = {key: float(copy.deepcopy(self.dig_dataset[key].attrs["start_frequency"])) for key in self.keys} #starting frequencies of scans
		self.fstop = {key: float(copy.deepcopy(self.dig_dataset[key].attrs["stop_frequency"])) for key in self.keys} #ending frequencies of scans
		self.fres = {key: float(copy.deepcopy(self.dig_dataset[key].attrs['frequency_resolution'])) for key in self.keys}
		self.Q = {key: float(copy.deepcopy(self.dig_dataset[key].attrs["Q"])) for key in self.keys} #quality factor during scan
		self.notes = {key: copy.deepcopy(self.dig_dataset[key].attrs["notes"]) for key in self.keys} #notes attached to scan
		self.errors = {key: copy.deepcopy(self.dig_dataset[key].attrs['errors']) for key in self.keys} #errors attached to scan
		self.int_times = {key: float(copy.deepcopy(self.dig_dataset[key].attrs['integration_time'])) for key in self.keys}
		self.bad_scans = {key: self.dig_dataset[key] for key in self.keys if self.dig_dataset[key].attrs['cut']} #Dict containing bad datasets
		
		add_input(self.dig_dataset,self.Tsys,'Tsys')
		if args!=None and self.clear_grand_spectra and 'grand_spectra_run1a' in self.h5py_file.keys():
			del self.h5py_file['grand_spectra_run1a']
		# derive necessary digitization structures??
		
		#Import axion frequencies to inject
		
		with open(axion_frequencies_filename, 'r') as file:
			self.axion_frequencies_to_inject = list(map(float, file.readlines()))
		
		# metadata (bad scans, etc.)

        # other definitions
		self.bad_scan_criteria = {'Tsys_bounds': (0.2,5.0),                							#Temperature bounds (lower, upper)
									'freq_bounds':(640.0,685.0),							#Frequency bounds (lower, upper). Bounds placed by run1a definition
									'Q_bounds':(10000,70000),                          				#Q bounds (lower, upper). Bounds retrieved from run1a_definitions.yaml
									'bad_timestamps': self.h5py_file["bad_timestamps_run1a"][...], #array-like
									'good_notes':("nan","filled_in_after_SAG"),
									'bad_logging': self.no_axion_log,
									'good_errors': ("", 
												'nan', 
												np.nan,
												' cannot get spectrum while acquiring',
												' cannot get spectrum while acquiring ca',
												' cannot start run, already acquiring',
												' cannot start run, already acquiring ca',
												'SetPowerupDefaultsPX4 in setup(): Phase'
												) 													#This tuple of errors is ok. If scan has any of these errors, analysis will still be performed
									 }
				 
		
		
		
		return None

    # defs for making decisions in the init

    # methods for calling executables

	def execute(self, timeit=False):
		try:
			bad_scans = self.collect_bad_scans()
			#Catch zero scan analysis
			if len(self.keys)==0 or not bad_scans:
				print("No scans to analyze")
				self.h5py_file.close()
				try:
					sys.exit(0)
				except SystemExit:
					os._exit(0)
			
			# set all calculations in motion			
			
			
			import oo_analysis.signals
			self.signal_dataset = oo_analysis.signals.generate(self)			
			
			import oo_analysis.analysis
			
			self.analysis_start=time.time()
			self.grand_spectra_group, self.candidates, ncut = oo_analysis.analysis.grand_spectra(self)
			self.analysis_stop=time.time()
			
			
			
			#import MCMC
			# perform MCMC analysis
			
			
			#import analytics
			


		except (KeyError, TypeError, SyntaxError) as error:
			self.h5py_file.close() #prevent corruption on break
			print("Execution failed with error: \n {0}".format(error))
			open(self.error_file, 'a+').write("\n\n"+ str(error))
			raise
		except KeyboardInterrupt:
			self.h5py_file.close()
			print("Interrupted")
			try:
				sys.exit(0)
			except SystemExit:
				os._exit(0)
		finally:
			#save analysis to disk, generate plots, and close out file
			string="Analysis of {0} scans took {1:0.3f} seconds. \tOf those scans, {2:d} were cut".format(len(self.keys),  self.analysis_stop-self.analysis_start, ncut)
			print(string)
			self.output()
			if self.make_plots:
				self.generate_plots()
			self.h5py_file.close() #Close the file, saving the changes.
			return None
			
	def output(self):
		# outputs data to local "./output/" directory
		import oo_analysis.data_management
		print("\nWriting output to:\n\t {0}".format(self.output_file))
		oo_analysis.data_management.write_out(self.grand_spectra_group,self.output_file)
		return None
		
	def collect_bad_scans(self):
		# runs through bad scan critereon and removes bad scans
		# collecting metadata for later analysis
		# may want to set up to run through only one or a set of conditions
		# should also try to make dynamics so that it only tries new conditions
		print("Collecting bad scans                        \r", end='')
		
		#Get bad scan criteria
		condition = self.bad_scan_criteria
		Q_bounds = condition["Q_bounds"]
		freq_bounds = condition["freq_bounds"]
		Tsys_bounds = condition["Tsys_bounds"]
		bad_timestamps = condition["bad_timestamps"]
		good_notes = condition["good_notes"]
		bad_logging = condition["bad_logging"]
		good_errors = condition["good_errors"]
		counter = 1
		N_iter = len(self.keys)
		bad_scans = []
		contains_good_scans=False
		for key in self.keys: # python 3.x
			try:
				cut = False
				cut_reason = ""
				
									
				#Get scan attributes
				timestamp = self.timestamp[key]
				Q = self.Q[key]
				fstart = self.fstart[key]
				fstop = self.fstop[key]
				Tsys = self.Tsys[key]
				notes = self.notes[key]
				errs = self.errors[key]

				#Compare scan attributes and bad scan criteria
				if key in self.bad_scans:
					cut = True
					cut_reason = "Scan marked as bad scan in bad_scans.dat"
				elif timestamp in bad_timestamps:
					cut=True
					cut_reason = "Scan not suitable for analysis (Bad timestamp)"
					self.bad_scans.append(key)
				elif not ((Q>Q_bounds[0])&(Q<Q_bounds[1])): #If Q lies outside Q bounds, cut scan
					cut=True
					cut_reason = "exceeding Q bound"
					self.bad_scans.append(key)
				elif ((fstart<freq_bounds[0])|(fstart>freq_bounds[1])) or ((fstop<freq_bounds[0])|(fstop>freq_bounds[1])):
					cut=True
					cut_reason = "scan outside 645 to 685 MHz range"
					self.bad_scans.append(key)
				elif ((Tsys>Tsys_bounds[1])|(Tsys<Tsys_bounds[0])):
					cut=True
					cut_reason = "exceeding system temperature bound"
					self.bad_scans.append(key)
				elif notes not in good_notes and not pd.isnull(notes):
					cut=True
					cut_reason = "Scan not suitable for analysis (See dataset notes)"
					self.bad_scans.append(key)
				elif key in bad_logging:
					cut=True
					cut_reason = "No associated axion log"
					self.bad_scans.append(key)
				elif errs not in good_errors and not pd.isnull(errs):
					cut = True
					cut_reason = "See errors"
					self.bad_scans.append(key)
				elif key in self.bad_scans:
					cut = True
					cut_reason = "Scan marked as bad scan in bad_scans.dat"
				if not cut:
					contains_good_scans=True
				self.dig_dataset[key].attrs["cut"] = cut
				self.dig_dataset[key].attrs["cut_reason"] = cut_reason
				bad_scans.append(cut)
			except (RuntimeError, KeyError) as error:
				print("\n\nError with scan {0}.".format(key))
				open(self.error_file, 'a+').write("\n\n"+ str(error))
				raise
			if int(counter/N_iter*100)-counter/N_iter*100<10**(-8):
				print("Collecting bad scans ({0} % complete) \r".format(int((counter/N_iter)*100)), end='')
			counter+=1
		return contains_good_scans

	def collect_meta_data(self):
		# collects summary data on each scan for decision-making and
		# later analysis
		return None
	
	def generate_plots(self):
		"""
		Method produces figures-of-merit from analysis.
		Current plots:
		
		coupling strength
		sensitivity to DM densities at coupling strength
		All plots in dissertation (coupling sensitivity, power sensitivity, DM density sensitivity)
		timestamp vs frequency (green when used, another color for thrown out scan for each reason (Q, Temperature, ...) Q wrt time
		placement of candidates (look up on elog)
		how many scans for certain frequency
		deltas as histogram (haystack for statistical things)
		plot of injected signal over many scans like what gray did
		different axion model injections (lentz, ITS, etc.), vary orbital modulations
		temperature wrt frequency and time
		"""
		
		from oo_analysis.figures.plotter import figures_class
		save_dir = os.getcwd() + "/oo_analysis/figures"
		kwargs = {"Tsys": list(self.Tsys.values()), "timestamp":list(self.timestamp.values()), "mode_frequencies": list(self.mode_frequencies.values()), "params": self.params}
		
		print("Generating plots")
		fig_cl = figures_class(save_dir)
		fig_cl.deltas(**kwargs)
		fig_cl.fit_significance(**kwargs)
		fig_cl.candidates()
		fig_cl.sensitivity_DM(**kwargs)
		fig_cl.sensitivity_power(**kwargs)
		fig_cl.sensitivity_coupling(**kwargs)
		fig_cl.SNR(**kwargs)
	
	def update_astronomical_tables(self):
		curr_time = dt.datetime.now()
		if "grand_spectra_run1a" in self.h5py_file.keys():
			last_change = self.h5py_file['grand_spectra_run1a']['last_change'][0].decode()
			if last_change!='':
				last_change_obj = dt.datetime.strptime(last_change, '%Y-%m-%d')
				if np.abs((curr_time-last_change_obj).total_seconds())/(60*60*24)>7.0: #if last execution was more than a week ago
					try:
						from astroplan import download_IERS_A
						print("\nUpdating astronomical tables\n")
						download_IERS_A()
						return None
					except urllib.error.URLError:
						print("Cannot download astronomical data. Check internet connection or manually update datatables. Defaulting to stored data")
	
	def garbage_collector(self):
		import gc
		gc.collect()
		