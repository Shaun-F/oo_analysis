"""
__init__.py: main init file for signal generation

Created by: Erik Lentz
Creation Date: 10/26/18
"""
from .modulation import modulation
from oo_analysis.toolbox.dict_to_object import dict_to_object


def generate(object):
    # get structures
	sigkeys = ["signal", "axion_mass", "nbins", "keys"]       # take down signal keys, like signal type and params
	sigparams = {key: getattr(object,key) for key in sigkeys}
	
	scankeys = ["start_scan", "end_scan", "dig_dataset", "timestamp", "fstart", "fstop", "Tsys"]     # take down keys for getting scan specific information, like scan number and central mass/frequency, timestep, etc
	scanparams_dict =  {key: getattr(object, key) for key in scankeys}
	scanparams = dict_to_object(**scanparams_dict)   # Change dictionary to object
	
	modkeys = ["pec_vel"]  # take down modulation keys, like which elements to include
	modparams = {key: getattr(object,key) for key in modkeys}
	
	# get callable signal
	
	# modulate signal for each scan
	signals = modulation(**modparams, **scanparams_dict, **sigparams).executor()
	
	return signals