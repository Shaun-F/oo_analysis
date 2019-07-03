"""
backsub_filters_lib.py: library for background subtraction operations

Created by Erik Lentz
Creation Date: 7/19/18
"""
import sys
sys.path.append("..")
def poly_fit(inputs,degree=6):
	import numpy
	signal = inputs
    # subtracts fitted polynomial from input arrays
    # retreive fitting kernel

    # execute multiprocessing routine for fitting

    # execute multiprocessing routine for subtraction
	domain = numpy.arange(len(inputs))
	polyfit = numpy.polyfit(domain, signal, degree)
	outputs = signal/(numpy.poly1d(polyfit)(domain))
	return outputs # in dictionary form

def SG(inputs,window=101,degree=4):
	from scipy.signal import savgol_filter as SGfilter
	signal = inputs
    # subtracts convolved polynomial from input arrays
    # fit using python packages

    # execute multiprocessing routine for subtraction
	
	to_remove = SGfilter(signal, window, degree)
	outputs = signal/to_remove
	return outputs # in dictionary form

def RCHPF(inputs,window=10,copies=3,dyn=False, **submeta):
	""" the only one we need right now
	import functinality from your RCHPF file
	"""
	from filters.RCHPF import DFT, IDFT, reciprocated_clone_hpf
	# subtracts large scale structure from input arrays
	

	sub_data = reciprocated_clone_hpf(inputs, copies, False, '', **submeta)
	# determine if dynamic and set flags

	
	########For distributed architecture.
	# retreive fitting kernel

	# execute multiprocessing routine for fitting

	# execute multiprocessing routine for subtraction


	return sub_data  #outputs # in dictionary form
