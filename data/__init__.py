"""
__init__.py: main init file for data

Created by: Erik Lentz
Creation Date: 10/26/18
"""
import sys
sys.path.append("../experiment")
sys.path.append("..")
from get_squid_dataset import get_squid_dataset
def input(parameters):
	# sets up start and stop parameters
	start = int(parameters['start_scan'])
	stop = int(parameters['end_scan'])
	# gets intput data
	import h5py
	import warnings
	warnings.simplefilter(action='ignore', category=FutureWarning)
	print("Loading hdf5 file and datasets")
	data_file = h5py.File(b"../data/raw/run1a_data.hdf5", "r+")
	
	if "bad_timestamps_run1a" not in data_file.keys():
		data_file.create_dataset(name="bad_timestamps_run1a", dtype="S10", data = [b"initval"], maxshape=(None,))
		
	dig_dataset = {}
	no_axion_log = []
	
	try:
		dig_dataset = {str(key): data_file['digitizer_log_run1a'][str(key)] for key in range(start, stop)}
		
	except (KeyError, MemoryError):
		raise
		
	for key in dig_dataset:
		if not 'alog_timestamp' in dig_dataset[key].attrs:
			no_axion_log.append(key)
	
	
	return dig_dataset, data_file, no_axion_log

def add_input(database,trait,trait_name):
	"""function takes trait in dictionary form and inputs attribute into
	 database stratified with same key
	 If dataset doesnt havent attribute <trait_name>, function automatically creates and populates it with <trait>
	"""
	for key in trait.keys():
		if trait_name in database[key].attrs:
			database[key].attrs[trait_name] = trait[key]
		elif trait_name not in database[key].attrs:
			database[key].attrs.create(trait_name, trait[key])
	return database
	
	"""
	for i in range(start,stop):
		try:
			inx = str(i)
			dig_dataset[inx] = data_file["digitizer_log_run1a"][inx]
		except KeyError:
			continue #No digitizer log found. check next scan number
		try:
			inx = str(i)
			axion_dataset[inx] = data_file["axion_log_run1a"][inx]
		except KeyError:
			no_axion_log.append(str(i))
			continue
	"""
