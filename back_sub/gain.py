"""
gain.py: holds routines for background subtraction code based on pure gain model

Created by: Erik Lentz
Creation Date: 10/26/18
"""
import filters.backsub_filters_lib as bslib


def execute(paramdict,datadict):
	deltas = {}
	# cyle though data
	for key,scan in datadict['dig_dataset'].items(): # or whatever the key is
		# get appropriate BS method (using RCHPF for placeholder)
		window = paramdict["filter_params"][0]
		copies = paramdict["filter_params"][1]
		scan_deltas = bslib.RCHPF(scan,window=window,copies=copies,dyn=False) #supply kwargs from paramdict
		deltas[key] = scan_deltas

	return deltas
