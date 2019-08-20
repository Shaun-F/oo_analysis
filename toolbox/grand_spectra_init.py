def initialization(chunk, scan_keys):
	#Initialize all datasets with empty containers
	import numpy
	datanames_zeroed=['sigma_w', 'optimal_weight_sum', 'noise_power', 'model_excess_sqrd', 'axion_fit_significance', 'power_deviation', 'nscans', 'SNR']
	stringnames=['scans', 'scans_in', 'scans_out']
	datanames_inf = ['axion_fit_uncertainty', 'sensitivity_power', 'sensitivity_coupling','axion_fit']
	
	
	init_data = numpy.arange(start = 644*10**6, stop = 682*10**6, step = 95.4)
	keys = chunk.keys()
	for i in datanames_zeroed:
		if i not in keys:
			chunk.create_dataset(name=str(i), dtype=numpy.float64, data=numpy.asarray([0]*len(init_data), dtype='byte'), maxshape=(None,))
		else:
			pass
	for i in datanames_inf:
		if i not in keys:
			chunk.create_dataset(name=str(i), dtype=numpy.float64, data=numpy.asarray([numpy.inf]*len(init_data)), maxshape=(None,))
		else:
			pass
	for i in stringnames:
		if i not in keys:
			chunk.create_dataset(name=str(i), dtype="S10", data=numpy.asarray([b""]*170000), maxshape=(None,))
		else:
			pass
	
	if 'last_change' not in keys:
		chunk.create_dataset(name='last_change', dtype="S10", data = numpy.asarray([b""]), maxshape=(None,))
	
	if 'axion_frequencies' not in keys:
		chunk.create_dataset(name='axion_frequencies', dtype = numpy.float64, data = init_data, maxshape = (None,))
	
	chunk.create_group('deltas')
	[chunk['deltas'].create_dataset(name=str(key), dtype=numpy.float64, shape = (256,), maxshape = (None,)) for key in scan_keys]
	