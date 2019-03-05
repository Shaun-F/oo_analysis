"""
core.py: holds  control routines for analysis

Created by: Erik Lentz
Creation Date: 6/1/18
"""
import sys
sys.path.append("../experiment")
sys.path.append("..")
from control.param_parser import parser
from experiment.get_squid_dataset import get_squid_dataset
from experiment.calc_sys_temp_offline import calc_sys_temp
import time
import argparse


############# Argument parsing
P = argparse.ArgumentParser(description="Main execution file for oo_analysis")
P.add_argument('-t', '--timeit', action='store', default=False, help='Argument specifies whether to time all the subprocesses of the analysis')
args = P.parse_args()

timeit = args.timeit
#############

# create main structure
class core_analysis():
	def __init__(self):
        # get parameters
		filename = "../job.param"  # from command line, nominally
		params = parser(filename)
        # set switches of analysis
        # default
        # from param file
		for key,value in params.items(): # for python 3.x
			setattr(self,key,value)
		
        # set data structures
		import data.__init__ 
		
		pulldata_start = time.time()
		self.dig_dataset, self.h5py_file, self.no_axion_log = data.__init__.input(params)
		self.keys = [i for i in self.dig_dataset.keys() if i not in self.no_axion_log] #Was originally dig_dataset,but some digitizer logs didnt have associated axion logs.
		pulldata_stop = time.time()
		
		self.output_file = "../output/grand_spectra.dat"
		print("Loading data successful. It took {0:0.3f} seconds. Beginning analysis of {1} scans".format((pulldata_stop-pulldata_start), len(self.dig_dataset)))
		
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
			
			self.collect_bad_scans()
			import signals
			signals_start = time.time()
			self.signal_dataset = signals.generate(self)
			signals_stop=time.time()
			
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
				self.meta_analysis.append("Total analysis time is {0:03f} seconds".format(analysis_stop - analysis_start))
				[print(x) for x in self.meta_analysis[1:]]
				print("###############################################################################\n\n")
				with open("../meta/analysis_statistics.txt", "w") as f:
					for x in self.meta_analysis[1:]:
						f.write(str(x) + "\n")
			self.output()
			self.h5py_file.close() #Close the file, saving the changes.
			
			string=" \n Analysis of {0} scans took {1:0.3f} seconds. \n Of those scans, {2:d} were cut. \n\n Writing output to {3}".format(len(self.keys),  analysis_stop-analysis_start, ncut, self.output_file)
			
			print(string)
			return None
		except (KeyError, TypeError) as error:
			self.h5py_file.close() #prevent corruption on break
			print(error)
			raise
			



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
			self.dig_dataset[key].attrs["cut"] = cut
			self.dig_dataset[key].attrs["cut_reason"] = cut_reason
			bad_scans_stop = time.time()
			bad_scans_timer.append(bad_scans_stop-bad_scans_start)
		
		#meta analysis
		if self.meta_analysis[0]:
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