"""
__init__.py: main init file for the background subtraction routines

Created by: Erik Lentz
Creation Date: 10/26/18
"""
# set up for background subtraction
import back_sub.gain as gain
import sys
sys.path.append("..")
from param_parser import parser
from toolbox.pulldata import pulldata

def BS(object):
	# get structures
	file_name = "../job.param"
	params = parser(file_name)
	attributes = dir(object)
	bskeys = ["filter", "filter_params", "parallel", "signal", "pec_vel"] # take down basic keys, like filter type and parameters
	bsparams = {key: getattr(object,key) for key in bskeys}
	datakeys = ["pull_data", "start_scan", "end_scan", "bad_scans", "temp", 'scans']# take down data keys
	if object.pull_data == 'false':
		pass
	elif object.pull_data == 'true':
		pulldata(object) #pull data from file and attach to object
	bsdata = {key: getattr(object,key) for key in datakeys}
	# run routine
	deltas = gain.execute(bsparams,bsdata)
	# set bs'ed data in object
	
	[setattr(object,key,deltas[key]) for key in deltas.keys()] # may want to nest in original dataset or not
	return object
