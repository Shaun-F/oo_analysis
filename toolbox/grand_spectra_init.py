def initialization(chunk):
	#Initialize all datasets with empty containers
	import numpy
	datanames_zeroed=['sigma_w', 'optimal_weight_sum', 'noise_power', 'model_excess_sqrd','axion_fit', 'axion_fit_significance', 'power_deviation', 'nscans']
	stringnames=['scans', 'scans_in', 'scans_out', 'last_change']
	datanames_inf = [ 'axion_fit_uncertainty', 'sensitivity_power', 'sensitivity_coupling','SNR']
	
	
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
	
	if 'axion_frequencies' not in keys:
		chunk.create_dataset(name='axion_frequencies', dtype = numpy.float64, data = init_data, maxshape = (None,))
	
	
	
	"""
	if 'scans' not in chunk.keys():
		chunk.create_dataset(name='scans', dtype='byte', data=numpy.asarray([]), maxshape=(None,)) #array of scan id strings
	if 'sigma_w' not in chunk.attrs:
		chunk.create_dataset(name='sigma_w', dtype=numpy.float64, data=numpy.asarray([]), maxshape=(None,))
	if 'optimal_weight_sum' not in chunk.attrs:
		chunk.create_dataset(name='optimal_weight_sum', dtype=numpy.float64, data=numpy.asarray([]), maxshape=(None,))
	if 'SNR' not in chunk.attrs:
		chunk.create_dataset(name='SNR', dtype=numpy.float64, data=numpy.asarray([]), maxshape=(None,))
	if 'noise_power' not in chunk.attrs:
		chunk.create_dataset(name='noise_power', dtype=numpy.float64, data=numpy.asarray([]), maxshape=(None,))
	if 'model_excess_sqrd' not in chunk.attrs:
		chunk.create_dataset(name='model_excess_sqrd', dtype=numpy.float64, data=numpy.asarray([]), maxshape=(None,))
	if 'axion_fit' not in chunk.attrs:
		chunk.create_dataset(name='axion_fit', dtype=numpy.float64, data=numpy.asarray([]), maxshape=(None,))
	if 'axion_fit_uncertainty' not in chunk.attrs:
		chunk.create_dataset(name
		chunk.attrs['axion_fit_uncertainty'] = []
	if 'sensitivity_power' not in chunk.attrs:
		chunk.attrs['sensitivity_power'] = []
	if 'sensitivity_coupling' not in chunk.attrs:
		chunk.attrs['sensitivity_coupling'] = []
	if 'axion_frequencies' not in chunk.attrs:
		chunk.attrs['axion_frequencies'] = []
	if 'power_deviation' not in chunk.attrs:
		chunk.attrs["power_deviation"] = []
	if 'nscans' not in chunk.attrs:
		chunk.attrs['nscans'] = []
	if 'scans_in' not in chunk.attrs:
		chunk.attrs['scans_in'] = np.asarray([], dtype='byte') #array of scan id strings in grand spectra
	if 'scans_out' not in chunk.attrs:
		chunk.attrs['scans_out'] = np.asarray([], dtype='byte') #array of scan id strings not in grand spectra
	if 'last_change' not in chunk.attrs:
		chunk.attrs['last_change'] = np.asarray([], dtype='byte') #array of timestamp strings
		
	"""