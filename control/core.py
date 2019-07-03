"""
core.py: holds  control routines for analysis

Created by: Erik Lentz and Modified by Shaun Fell
Creation Date: 6/1/18
"""
import sys; import os
sys.path.append("../experiment")
sys.path.append("..")
from control.param_parser import parser
from experiment.get_squid_dataset import get_squid_dataset
from experiment.calc_sys_temp_offline import calc_sys_temp
import time; import datetime
import argparse
import copy 
from toolbox.plot_dataset import plotter
from toolbox.freq_to_mass import freq_to_mass
from analysis.synthetic_injection import axion_injector

############# Argument parsing
P = argparse.ArgumentParser(description="Main execution file for oo_analysis")
P.add_argument('-t', '--timeit', action='store_true', default=False, help='Argument specifies whether to time all the subprocesses of the analysis')
P.add_argument('-cgs', '--clear_grand_spectra', action='store_true', default=False, help="Argument specifies whether to delete the grand spectra and start from scratch. Default is False")
P.add_argument('--start_scan', action='store', default = '', help="Argument specifies the starting scan number of the analysis. If not specified, starting number specified by job.param")
P.add_argument('--end_scan', action='store', default = '', help="Argument specifies the ending scan number of the analysis. If not specified, ending number specified by job.param")
P.add_argument('-p', '--make_plots', action='store', default=False, help="Argument specifies whether to generate plots or not after analysis is completed")

args = P.parse_args()

timeit = args.timeit
reset = args.clear_grand_spectra
start_scan = args.start_scan
end_scan = args.end_scan
make_plots = args.make_plots
#############

# create main structure
class core_analysis():
	def __init__(self, **kwargs):
        # get parameters
		filename = "../job.param"  # from command line, nominally
		params = parser(filename)
        # set switches of analysis
        # default
        # from param file

		if start_scan!='':
			params['start_scan']=start_scan
		if end_scan!='':
			params['end_scan']=end_scan
		
		for arg, val in kwargs.items():
			if arg in params.keys():
				params[str(arg)] = val
				
		for key,value in params.items(): # for python 3.x
			setattr(self,key,value)
			
			#find bad scans saved to file
		self.bad_scans_file = open(str(params['bad_scans_file']))
		self.bad_scans = self.bad_scans_file.read().splitlines()
		params['bad_scans'] = self.bad_scans	
			
        # set data structures
		import data.__init__ 
		pulldata_start = time.time()
		self.dig_dataset, self.h5py_file, self.no_axion_log, self.partitioned = data.__init__.input(params)
		self.keys = [copy.deepcopy(i) for i in self.dig_dataset.keys() if i not in self.no_axion_log] #Was originally dig_dataset,but some digitizer logs didnt have associated axion logs.
		pulldata_stop = time.time()
		
		
		self.output_file = "../output/grand_spectra.dat"
		print("Loading data successful. It took {0:0.3f} seconds. Beginning analysis of {1} scans".format((pulldata_stop-pulldata_start), len(self.keys)))
		if self.partitioned:
			print("Note: Dataset was partitioned")
		if reset and 'grand_spectra_run1a' in self.h5py_file.keys():
			del self.h5py_file['grand_spectra_run1a']	
        
		
		# derive necessary experiment data structures (put into dig_dataset)
		#Populate parameter dict's with dataset attributes.
		self.Tsys = {key: copy.deepcopy(self.dig_dataset[key].attrs["squid_temperature"]) for key in self.keys} #temperature of system during scan
		self.timestamp = {key: copy.deepcopy(self.dig_dataset[key].attrs["alog_timestamp"]) for key in self.keys} #timestamp of scan
		self.mode_frequencies = {key: copy.deepcopy(self.dig_dataset[key].attrs['mode_frequency']) for key in self.keys} #mode frequency of scan
		self.axion_mass = {key: freq_to_mass(copy.deepcopy(self.mode_frequencies[key])*10**6) for key in self.keys} #axion mass in eV with corresponding frequency equal to mode frequency
		self.fstart = {key: float(copy.deepcopy(self.dig_dataset[key].attrs["start_frequency"])) for key in self.keys} #starting frequencies of scans
		self.fstop = {key: float(copy.deepcopy(self.dig_dataset[key].attrs["stop_frequency"])) for key in self.keys} #ending frequencies of scans
		self.Q = {key: float(copy.deepcopy(self.dig_dataset[key].attrs["Q"])) for key in self.keys} #quality factor during scan
		self.notes = {key: copy.deepcopy(self.dig_dataset[key].attrs["notes"]) for key in self.keys} #notes attached to scan

		data.add_input(self.dig_dataset,self.Tsys,'Tsys')
		
		# Inject synthetic axion signals into datasets
		axion_injector(self)
		
		
		# derive necessary digitization structures??
		
		# metadata (bad scans, etc.)

        # other definitions
		self.bad_scan_criteria = {'Tsys': 5.0,                                                               #Place holder value. Update bound?
         'timestamps': self.h5py_file["bad_timestamps_run1a"][...], #array-like
         'freq_low':644.9, 'freq_top':680.1,
         'notes_neq':("nan","filled_in_after_SAG"),
		 'Q':10**6,
		 'bad_logging': self.no_axion_log
		 }
		 
		
		
		
		return None

    # defs for making decisions in the init

    # methods for calling executables

	def execute(self, timeit=False):
		try:
			#Catch zero scan analysis
			if len(self.keys)==0:
				print("No scans to analysis")
				try:
					sys.exit(0)
				except SystemExit:
					os._exit(0)
			# sets all calculations in motion
			self.meta_analysis = [timeit]
			"""
			import back_sub.__init__
			self = back_sub.__init__.BS(self)
			self.bad_scan_criteria['background'] = 'background condition'
			"""
			self.collect_bad_scans()
			print("generating signals")
			import signals
			self.signals_start = time.time()
			self.signal_dataset = signals.generate(self)
			self.signals_stop=time.time()
			
			print("signal generating complete. Beginning Analysis")
			import analysis
			self.analysis_start=time.time()
			self.grand_spectra_group, ncut = analysis.grand_spectra(self)
			self.analysis_stop=time.time()
			
			
			#import MCMC
			# perform MCMC analysis
			
			
			#import analytics
			
			
			# generate plots

		except (KeyError, TypeError, SyntaxError) as error:
			self.h5py_file.close() #prevent corruption on break
			print("Execution failed with error: \n {0}".format(error))
			open('../meta/error_log', 'a+').write(str(time.time())+ "\n\n"+ str(error))
			raise
		except KeyboardInterrupt:
			self.h5py_file.close()
			print("Interrupted")
			try:
				sys.exit(0)
			except SystemExit:
				os._exit(0)
		finally:
			#save analysis to disk and close out file
			string=" \n Analysis of {0} scans took {1:0.3f} seconds. \n Of those scans, {2:d} were cut".format(len(self.keys),  self.analysis_stop-self.analysis_start, ncut)
			print(string)
			self.collect_bad_scans()
			self.collect_meta_data()
			self.output()
			#self.generate_plots()
			self.h5py_file.close() #Close the file, saving the changes.
			return None
			
	def output(self):
		# outputs data to local "./output/" directory
		import data_management
		print("\n\nWriting output to {0}".format(self.output_file))
		data_management.write_out(self.grand_spectra_group,self.output_file)
		return None
		
	def collect_bad_scans(self):
		# runs through bad scan critereon and removes bad scans
		# collecting metadata for later analysis
		# may want to set up to run through only one or a set of conditions
		# should also try to make dynamics so that it only tries new conditions
		bad_scans_timer = []
		for key in self.keys: # python 3.x
			try:
				if key not in self.bad_scans:
					bad_scans_start = time.time()
					cut = False
					cut_reason = ""
					condition = self.bad_scan_criteria
					if self.timestamp[key] in condition["timestamps"]:
						cut=True
						cut_reason = "Scan not suitable for analysis (Bad timestamp)"
						self.bad_scans.append(key)
					if self.Q[key]>condition["Q"]:
						cut=True
						cut_reason = "exceeding Q bound"
						self.bad_scans.append(key)
					elif self.fstart[key]<condition["freq_low"] or self.fstop[key]<condition["freq_low"] or self.fstart[key]>condition["freq_top"] or self.fstop[key]>condition["freq_top"]:
						cut=True
						cut_reason = "scan outside 645 to 685 MHz range"
						self.bad_scans.append(key)
					elif self.Tsys[key]>condition["Tsys"]:
						cut=True
						cut_reason = "System temperature over bound"
						self.bad_scans.append(key)
					elif self.notes[key] in condition["notes_neq"]:
						cut=True
						cut_reason = "Scan not suitable for analysis (See dataset notes)"
						self.bad_scans.append(key)
					elif key in condition['bad_logging']:
						cut=True
						cut_reason = "No associated axion log"
						self.bad_scans.append(key)
					self.dig_dataset[key].attrs["cut"] = cut
					self.dig_dataset[key].attrs["cut_reason"] = cut_reason
					bad_scans_stop = time.time()
					bad_scans_timer.append(bad_scans_stop-bad_scans_start)
			except (RuntimeError, KeyError) as error:
				print("\n\nError with scan {0}.".format(key))
				open('../meta/error_log', 'a+').write(str(time.time())+ "\n\n"+ str(error))
				raise
		#meta analysis
		if self.meta_analysis[0]:
			if len(bad_scans_timer)==0:
				average_time = 0
			else:
				average_time = sum(bad_scans_timer)/len(bad_scans_timer)
			self.meta_analysis.append("Collecting bad scans took {0:0.3f} seconds".format(average_time))
		return None

	def collect_meta_data(self):
		# collects summary data on each scan for decision-making and
		# later analysis
		if self.meta_analysis[0]:
			print("\n\n################################ Meta-analysis ################################")
			self.meta_analysis.append("\n\nTotal signal generation time is {0:03f} seconds".format(self.signals_stop - self.signals_start))
			self.meta_analysis.append("Total analysis time of {0} scans is {1:03f} seconds".format(len(self.keys), self.analysis_stop - self.analysis_start))
			[print(x) for x in self.meta_analysis[1:]] #zeroth index contains parameters of meta-analysis
			print("#################################### End of Meta-analysis###########################\n\n")
			date = datetime.datetime.now()
			filename = "../meta/analysis_statistics" + "(" + str(date.year) + "-" + str(date.month) + "-" + str(date.day) + "_" + str(date.hour) + "-" + str(date.minute) + "-" + str(date.second) + ")(" + str(len(self.keys)) + "_scans).txt"
			with open(filename, "w") as f:
				for x in self.meta_analysis[1:]:
					f.write(str(x) + "\n")

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
		
		import figures.plotter; 
		from figures.plotter import figures_class
		save_dir = "C:/users/drums/documents/coding software/python/scripts/New-Analysis-Scheme/oo_analysis/figures"
		kwargs = {"Tsys": list(self.Tsys.values()), "timestamp":list(self.timestamp.values()), "mode_frequencies": list(self.mode_frequencies.values())}
		
		fig_cl = figures_class(save_dir)
		fig_cl.temp_v_freq(**kwargs)
		fig_cl.temp_v_time(**kwargs)
		fig_cl.deltas(**kwargs)
		fig_cl.candidates()
		fig_cl.sensitivity_DM()
		fig_cl.sensitivity_power()
		fig_cl.sensitivity_coupling()
		fig_cl.SNR()
		
		
		
if __name__ in '__main__':
	total_analysis_start = time.time()
	x = core_analysis()
	x.execute(timeit)
	n=1
	
	while x.partitioned:
		n+=1
		reset=False
		print("running next partition. (loop number: {0}) (partition [{1}, {2}])".format(n, max(map(int,x.keys))+1, x.end_scan))
		params = {'start_scan': max(map(int,x.keys))+1, 'end_scan': x.end_scan}
		x = core_analysis(**params)
		x.execute(timeit)
	if make_plots:
		x.generate_plots()
		
	total_analysis_stop = time.time()
	print("\n\n Entire analysis over all {0} partitions took {1:0.3f} seconds".format(n, total_analysis_stop-total_analysis_start))
	
		
		
		
		
		
		
		
		
		
		