"""
__init__.py: main init file for the background subtraction routines

Created by: Erik Lentz
Creation Date: 10/26/18
"""
# set up for background subtraction
import back_sub.gain as gain
import sys;sys.path.append("..")
from param_parser import parser
from toolbox.pulldata import pulldata

def BS(object):
	# get structures
		
	bskeys = ["filter", "filter_params", "parallel", "signal", "pec_vel"] # take down basic keys, like filter type and parameters
	bsparams = {key: getattr(object,key) for key in bskeys} #Pull keys from object class
	datakeys = ["pull_data", "start_scan", "end_scan", "bad_scans", "temp", 'dig_dataset']# take down data keys
	bsdata = {key: getattr(object,key) for key in datakeys} #Pull keys from object class
	# run routine
	deltas = gain.execute(bsparams,bsdata, object.meta_analysis)
	# set bs'ed data in object
	
	[setattr(object,key,deltas[key]) for key in deltas.keys()] # may want to nest in original dataset or not
	return object
