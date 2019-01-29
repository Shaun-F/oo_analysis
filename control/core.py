"""
core.py: holds  control routines for analysis

Created by: Erik Lentz
Creation Date: 6/1/18
"""
from param_parser import parser
import sys
sys.path.append("../experiment")
sys.path.append("..")
from get_squid_dataset import get_squid_dataset
from calc_sys_temp_offline import calc_sys_temp
import time

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
		self.dig_dataset, self.axion_dataset, self.squid_dataset, self.h5py_file = data.__init__.input(params)
		self.keys = self.dig_dataset.keys()
		
        # derive necessary experiment data structures (put into axion_dataset)
		#Populate parameter dict's with dataset attributes.
		self.Tsys = {key: calc_sys_temp(self.dig_dataset[key], self.axion_dataset[key], self.squid_dataset[key])["system temperature"] for key in self.keys}
		self.timestamp = {key: self.axion_dataset[key].attrs["timestamp"] for key in self.keys}
		self.fstart = {key: float(self.dig_dataset[key].attrs["start_frequency"]) for key in self.keys}
		self.fstop = {key: float(self.dig_dataset[key].attrs["stop_frequency"]) for key in self.keys}
		self.Q = {key: float(self.axion_dataset[key].attrs["Q"]) for key in self.keys}
		self.notes = {key: self.axion_dataset[key].attrs["notes"] for key in self.keys}
		self.afreq = self.axion_mass/(4.13566766*10**(-15)) #attach axion freq to object
			
		data.add_input(self.axion_dataset,self.Tsys,'Tsys')
		# derive necessary digitization structures??
		
		# metadata (bad scans, etc.)

        # other definitions
		self.bad_scan_criteria = {'Tsys': 5.0,                                                               #Place holder value. Update bound?
         'timestamps': self.h5py_file["bad_timestamps_run1a"][...], #array-like
         'freq_low':644.9, 'freq_top':680.1,
         'notes_neq':("nan","filled_in_after_SAG"),
		 'Q':10**6}
		return None

    # defs for making decisions in the init

    # methods for calling executables

	def execute(self):
		
		# sets all calculations in motion
		self.collect_bad_scans()
		import back_sub.__init__
		self = back_sub.__init__.BS(self)
		self.bad_scan_criteria['background'] = 'background condition'
		self.collect_bad_scans()
		import signals
		signals_start = time.time()
		self.signal_dataset = signals.generate(self)
		signals_stop=time.time()
		import analysis
		analysis_start=time.time()
		self.analysis_dataset = analysis.grand_spectra(self)
		analysis_stop=time.time()
		#import MCMC
		# perform MCMC analysis
		#import analytics
		# generate plots
		self.output()
		self.h5py_file.close() #Close the file, saving the changes.
		string="Signal generation took {0:0.3f} seconds. \n Analysis took {1:0.3f} seconds".format(signals_stop-signals_start, analysis_stop-analysis_start)
		print(string)
		return None



	def output(self):
		# outputs data to local "./output/" directory
		import data_management
		data_management.write_out(self.analysis_dataset,"../output/grand_spectra.dat")
		return None
		
	def collect_bad_scans(self):
		# runs through bad scan critereon and removes bad scans
		# collecting metadata for later analysis
		# may want to set up to run through only one or a set of conditions
		# should also try to make dynamics so that it only tries new conditions
		for key,value in self.dig_dataset.items(): # python 3.x
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
			self.axion_dataset[key].attrs["cut"] = cut
			self.axion_dataset[key].attrs["cut_reason"] = cut_reason
			self.dig_dataset[key].attrs["cut"] = cut
			self.dig_dataset[key].attrs["cut_reason"] = cut_reason
		return None

	def collect_meta_data(self):
		# collects summary data on each scan for decision-making and
		# later analysis

		return None

x = core_analysis()
x.execute()