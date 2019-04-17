"""
__init__.py: main init file for signal generation

Created by: Erik Lentz
Creation Date: 10/26/18
"""
import sys
sys.path.append("..")
import signals.signal_lib as slib
import signals.modulation as mod
from toolbox.dict_to_object import dict_to_object
import time


def generate(object):
    # get structures
	sigkeys = ["signal", "axion_mass", "nbins", "keys"]       # take down signal keys, like signal type and params
	sigparams = {key: getattr(object,key) for key in sigkeys}
	
	scankeys = ["start_scan", "end_scan", "dig_dataset", "timestamp", "fstart", "fstop", "Tsys"]     # take down keys for getting scan specific information, like scan number and central mass/frequency, timestep, etc
	scanparams_dict =  {key: getattr(object, key) for key in scankeys}
	scanparams = dict_to_object(**scanparams_dict)   # Change dictionary to object
	
	modkeys = ["pec_vel"]       # take down modulation keys, like which elements to include
	modparams = {key: getattr(object,key) for key in modkeys}
	
	# get callable signal
	signal_gen_start = time.time()
	signal = slib.signal(sigparams)
	signal_gen_stop = time.time()
	
	#meta analysis
	meta = {"meta_analysis":object.meta_analysis}
	
	# modulate signal for each scan
	signal_mod_start = time.time()
	signals = mod.modulation(**modparams, **scanparams_dict, **sigparams, **meta).executor()
	signal_mod_stop = time.time()
	
	"""
	if object.meta_analysis[0]:
		object.meta_analysis.append("Signal generation took {0:03f} seconds".format(signal_gen_stop - signal_gen_start))
		object.meta_analysis.append("Signal modulation took {0:03f} seconds".format(signal_mod_stop - signal_gen_start))
	"""
	return signals