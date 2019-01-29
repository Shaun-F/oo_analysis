"""
backsub_filters_lib.py: library for background subtraction operations

Created by Erik Lentz
Creation Date: 7/19/18
"""
import sys
sys.path.append("..")
def poly_fit(inputs,degree=6):
    # subtracts fitted polynomial from input arrays
    # retreive fitting kernel

    # execute multiprocessing routine for fitting

    # execute multiprocessing routine for subtraction

    return outputs # in dictionary form

def SG(inputs,window=101,degree=4):
    # subtracts convolved polynomial from input arrays
    # fit using python packages

    # execute multiprocessing routine for subtraction

    return outputs # in dictionary form

def RCHPF(inputs,window=10,copies=3,dyn=False):
	""" the only one we need right now
	import functinality from your RCHPF file
	"""
	from filters.RCHPF import DFT, IDFT, reciprocated_clone_hpf
	# subtracts large scale structure from input arrays
	sub_data = reciprocated_clone_hpf(data=inputs, npairs=copies)
	# determine if dynamic and set flags

	
	########For distributed architecture.
	# retreive fitting kernel

	# execute multiprocessing routine for fitting

	# execute multiprocessing routine for subtraction


	return sub_data  #outputs # in dictionary form
