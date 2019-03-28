"""
core.py: holds  control routines for analysis

Created by: Erik Lentz
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

############# Argument parsing
P = argparse.ArgumentParser(description="Main execution file for oo_analysis")
P.add_argument('-t', '--timeit', action='store', default=False, help='Argument specifies whether to time all the subprocesses of the analysis')
P.add_argument('-cgs', '--clear_grand_spectra', action='store', default=False, help="Argument specifies whether to delete the grand spectra and start from scratch. Default is False")
P.add_argument('--start_scan', action='store', default = '', help="Argument specifies the starting scan number of the analysis. If not specified, starting number specified by job.param")
P.add_argument('--end_scan', action='store', default = '', help="Argument specifies the ending scan number of the analysis. If not specified, ending number specified by job.param")
args = P.parse_args()

timeit = args.timeit
reset = args.clear_grand_spectra
start_scan = args.start_scan
end_scan = args.end_scan
#############

# create main structure
class core_analysis():
	def __init__(self, **kwargs):
        # get parameters
		filename = "../job.param"  # from command line, nominally
		params = parser(filename)
		for arg, val in kwargs.items():
			if arg in params.keys():
				params[str(arg)] = val
				
        # set switches of analysis
        # default
        # from param file
		for key,value in params.items(): # for python 3.x
			setattr(self,key,value)
		
        # set data structures
		import data.__init__ 
		
		if start_scan!='':
			params['start_scan']=start_scan
		if end_scan!='':
			params['end_scan']=end_scan
		pulldata_start = time.time()
		self.dig_dataset, self.h5py_file, self.no_axion_log, self.paritioned = data.__init__.input(params)
		self.keys = [i for i in self.dig_dataset.keys() if i not in self.no_axion_log] #Was originally dig_dataset,but some digitizer logs didnt have associated axion logs.
		pulldata_stop = time.time()
		
		self.output_file = "../output/grand_spectra.dat"
		print("Loading data successful. It took {0:0.3f} seconds. Beginning analysis of {1} scans".format((pulldata_stop-pulldata_start), len(self.keys)))
		if self.paritioned:
			print("Note: Dataset was paritioned")
		if reset and 'grand_spectra_run1a' in self.h5py_file.keys():
			del self.h5py_file['grand_spectra_run1a']	
		
        # derive necessary experiment data structures (put into dig_dataset)
		#Populate parameter dict's with dataset attributes.
		self.Tsys = {key: self.dig_dataset[key].attrs["squid_temperature"] for key in self.keys}
		self.timestamp = {key: self.dig_dataset[key].attrs["alog_timestamp"] for key in self.keys}
		self.fstart = {key: float(self.dig_dataset[key].attrs["start_frequency"]) for key in self.keys}
		self.fstop = {key: float(self.dig_dataset[key].attrs["stop_frequency"]) for key in self.keys}
		self.Q = {key: float(self.dig_dataset[key].attrs["Q"]) for key in self.keys}
		self.notes = {key: self.dig_dataset[key].attrs["notes"] for key in self.keys}
		#self.afreq = self.axion_mass/(4.13566766*10**(-15)) #attach axion freq to object

		data.add_input(self.dig_dataset,self.Tsys,'Tsys')
		# derive necessary digitization structures??
		
		# metadata (bad scans, etc.)

        # other definitions
		self.bad_scan_criteria = {'Tsys': 5.0,                                                               #Place holder value. Update bound?
         'timestamps': self.h5py_file["bad_timestamps_run1a"][...], #array-like
         'freq_low':644.9, 'freq_top':680.1,
         'notes_neq':("nan","filled_in_after_SAG"),
		 'Q':10**6,
		 'bad_logging': self.no_axion_log}
		return None

    # defs for making decisions in the init

    # methods for calling executables

	def execute(self, timeit=False):
		try:
			# sets all calculations in motion
			self.meta_analysis = [timeit]
			self.collect_bad_scans()
			
			"""
			import back_sub.__init__
			self = back_sub.__init__.BS(self)
			self.bad_scan_criteria['background'] = 'background condition'
			"""
			print("generating signals")
			self.collect_bad_scans()
			import signals
			signals_start = time.time()
			self.signal_dataset = signals.generate(self)
			signals_stop=time.time()
			
			print("signal generating complete. Beginning Analysis")
			import analysis
			analysis_start=time.time()
			self.analysis_dataset, ncut = analysis.grand_spectra(self)
			analysis_stop=time.time()
			
			#import MCMC
			# perform MCMC analysis
			#import analytics
			# generate plots
			if self.meta_analysis[0]:
				print("\n\n################################ Meta-analysis ################################")
				self.meta_analysis.append("\n\nTotal signal generation time is {0:03f} seconds".format(signals_stop - signals_start))
				self.meta_analysis.append("Total analysis time of {0} scans is {1:03f} seconds".format(len(self.keys), analysis_stop - analysis_start))
				[print(x) for x in self.meta_analysis[1:]]
				print("#################################### End of Meta-analysis###########################\n\n")
				date = datetime.datetime.now()
				filename = "../meta/analysis_statistics" + "(" + str(date.year) + "-" + str(date.month) + "-" + str(date.day) + "_" + str(date.hour) + "-" + str(date.minute) + "-" + str(date.second) + ")(" + str(len(self.keys)) + "_scans).txt"
				with open(filename, "w") as f:
					for x in self.meta_analysis[1:]:
						f.write(str(x) + "\n")
			self.output()
			self.h5py_file.close() #Close the file, saving the changes.
			
			string=" \n Analysis of {0} scans took {1:0.3f} seconds. \n Of those scans, {2:d} were cut. \n\n Writing output to {3}".format(len(self.keys),  analysis_stop-analysis_start, ncut, self.output_file)
			
			print(string)
			return None
		except (KeyError, TypeError, SyntaxError) as error:
			self.h5py_file.close() #prevent corruption on break
			print("Execution failed with error: \n {0}".format(error))
			raise
		except KeyboardInterrupt:
			self.h5py_file.close()
			print("Interrupted")
			try:
				sys.exit(0)
			except SystemExit:
				os._exit(0)
			



	def output(self):
		# outputs data to local "./output/" directory
		import data_management
		data_management.write_out(self.analysis_dataset,self.output_file)
		return None
		
	def collect_bad_scans(self):
		# runs through bad scan critereon and removes bad scans
		# collecting metadata for later analysis
		# may want to set up to run through only one or a set of conditions
		# should also try to make dynamics so that it only tries new conditions
		bad_scans_timer = []
		for key in self.keys: # python 3.x
			bad_scans_start = time.time()
			cut = False
			cut_reason = ""
			condition = self.bad_scan_criteria
			if self.timestamp[key] in condition["timestamps"]:
				cut=True
				cut_reason = "Scan not suitable for analysis (Bad timestamp)"
			if self.Q[key]>condition["Q"]:
				cut=True
				cut_reason = "exceeding Q bound"
			elif self.fstart[key]<condition["freq_low"] or self.fstop[key]<condition["freq_low"] or self.fstart[key]>condition["freq_top"] or self.fstop[key]>condition["freq_top"]:
				cut=True
				cut_reason = "scan outside 645 to 685 MHz range"
			elif self.Tsys[key]>condition["Tsys"]:
				cut=True
				cut_reason = "System temperature over bound"
			elif self.notes[key] in condition["notes_neq"]:
				cut=True
				cut_reason = "Scan not suitable for analysis (See dataset notes)"		
			elif key in condition['bad_logging']:
				cut=True
				cut_reason = "No associated axion log"
			self.dig_dataset[key].attrs["cut"] = cut
			self.dig_dataset[key].attrs["cut_reason"] = cut_reason
			bad_scans_stop = time.time()
			bad_scans_timer.append(bad_scans_stop-bad_scans_start)
		
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

		return None

		
if __name__ in '__main__':

	x = core_analysis()
	x.execute(timeit)
	while True:
		if x.paritioned:
			params = {'start_scan': max(x.keys), 'end_scan': x.end_scan}
			x = core_analysis(**params)
			x.execute(timeit)
			print("running next parition")
		if not x.paritioned:
			break
		
		
		
		
		
		
		
		
		
		
		