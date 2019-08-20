"""
gain.py: holds routines for background subtraction code based on pure gain model

Created by: Erik Lentz
Creation Date: 10/26/18
"""
import filters.backsub_filters_lib as bslib


def execute(paramdict,datadict, *meta):
	deltas = {}
	# cyle though data
	
			

	for key,scan in datadict['dig_dataset'].items(): # or whatever the key is
		# get appropriate BS method (using RCHPF for placeholder)
		window = paramdict["filter_params"][0]
		copies = paramdict["filter_params"][1]
		scan_deltas = bslib.RCHPF(scan[...],window,copies,False, **submeta) #supply kwargs from paramdict
		deltas[key] = scan_deltas["filtereddata"]
		submeta = scan_deltas['meta']
		
	#meta analysis
	if meta[0]:
		meta = meta[0] #remove surrounding tuple from *meta
		avg_reflecting_time = sum(submeta['reflecting_time'])/len(submeta['reflecting_time'])
		avg_calculating_tophat_size = sum(submeta['calculating_tophat_size'])/len(submeta['calculating_tophat_size'])
		avg_reciprocating_array = sum(submeta['reciprocating_array'])/len(submeta['reciprocating_array'])
		avg_generating_large_array = sum(submeta['generating_large_array'])/len(submeta['generating_large_array'])
		avg_generating_tophat = sum(submeta['generating_tophat'])/len(submeta['generating_tophat'])
		avg_fft_and_highpassfilter = sum(submeta['fft_and_highpassfilter'])/len(submeta['fft_and_highpassfilter'])
		avg_ifft = sum(submeta['ifft'])/len(submeta['ifft'])
		avg_picking_og_signal = sum(submeta['picking_og_signal'])/len(submeta['picking_og_signal'])
		avg_dividing_structure = sum(submeta['dividing_structure'])/len(submeta['dividing_structure'])
		meta.append("average time to reflect scan is {0:03f}".format(avg_reflecting_time))
		meta.append("average time to calculate tophat size is {0:03f}".format(avg_calculating_tophat_size))
		meta.append("average time to calculated reciprocated signal is {0:03f}".format(avg_reciprocating_array))
		meta.append("average time to generate large array of clones is {0:03f}".format(avg_generating_large_array))
		meta.append("average time to generate tophat is {0:03f}".format(avg_generating_tophat))
		meta.append("average time to fft and multiply by tophat is {0:03f}".format(avg_fft_and_highpassfilter))
		meta.append("average time to inverse fft is {0:03f}".format(avg_ifft))
		meta.append("average time to extract original signal is {0:03f}".format(avg_picking_og_signal))
		meta.append("average time to divide out structure is {0:03f}".format(avg_dividing_structure))
		

	return deltas
