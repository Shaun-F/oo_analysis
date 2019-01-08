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


def generate(object):
    # get structures
	sigkeys = ["signal", "axion_mass", "nbins"]       # take down signal keys, like signal type and params
	sigparams = {key: getattr(object,key) for key in sigkeys}
	
	scankeys = ["start_scan", "end_scan", "dig_dataset", "timestamp", "fstart", "fstop", "Tsys"]     # take down keys for getting scan specific information, like scan number and central mass/frequency, timestep, etc
	scanparams_dict =  {key: getattr(object, key) for key in scankeys}
	scanparams = dict_to_object(**scanparams_dict)   # Change dictionary to object
	
	modkeys = ["pec_vel"]       # take down modulation keys, like which elements to include
	modparams = {key: getattr(object,key) for key in modkeys}
	
	# get callable signal
	signal = slib.signal(sigparams)
	# modulate signal for each scan
	signals = mod.modulation(**modparams, **scanparams_dict, **sigparams).executor()
	"""
	for key,sparams in scanparams.__dict__.items(): # python 3.x
		signals[key] = mod.modulatedsignal(signal,scanparams,modparams)
	return signals
"""