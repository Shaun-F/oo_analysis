

def pulldata(core_analysis_object):
	"""
	Function takes core_analysis class, pulls data from .dig_dataset, and sets the new data as new attribute. 
	"""
	start = core_analysis_object.start_scan
	stop = core_analysis_object.end_scan
	dig_set = core_analysis_object.dig_dataset
	data = {}
	for x in range(start, stop):
		iter = str(x)
		data[iter] = dig_set[iter]
	
	setattr(core_analysis_object, 'scans', data)
	