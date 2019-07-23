"""
__init__.py: main init file for data

Created by: Erik Lentz
Creation Date: 10/26/18
"""
import time
import os
def input(parameters):
	# sets up start and stop parameters
	start = int(parameters['start_scan'])
	stop = int(parameters['end_scan'])
	# gets intput data
	import h5py
	import warnings
	warnings.simplefilter(action='ignore', category=FutureWarning)
	print("Loading hdf5 file and datasets")
	raw_data_filename = os.getcwd() + "/oo_analysis/data/raw/run1a_data.hdf5"
	data_file = h5py.File(raw_data_filename.encode(), "r+")
	
	if "bad_timestamps_run1a" not in data_file.keys():
		data_file.create_dataset(name="bad_timestamps_run1a", dtype="S10", data = [b"initval"], maxshape=(None,))
		
	dig_dataset = {}
	no_axion_log = []
	counter = 0
	partitioned = False
	
	try:
		for key in range(start, stop) :
			try:
				if key not in parameters['bad_scans']:
					dataset_toadd = data_file['digitizer_log_run1a'][str(key)]
					if 'alog_timestamp' in dataset_toadd.attrs:
						dig_dataset[str(key)] = dataset_toadd
						counter += 1
					if not 'alog_timestamp' in dataset_toadd.attrs:
						no_axion_log.append(str(key))
			except KeyError:
				pass
			if counter>=parameters['partition']:
				partitioned = True
				break

	except Exception as error:
		data_file.close()
		open('../meta/error_log', 'a+').write(str(time.time())+ "\n\n"+ str(error))
		raise
	return dig_dataset, data_file, no_axion_log, partitioned

def add_input(database,trait,trait_name):
	"""function takes trait in dictionary form and inputs attribute into
	 database stratified with same key
	 If dataset doesnt havent attribute <trait_name>, function automatically creates and populates it with <trait>
	"""
	for key in trait.keys():
		try:
		
			if trait_name in database[key].attrs:
				database[key].attrs[trait_name] = trait[key]
			elif trait_name not in database[key].attrs:
				database[key].attrs.create(trait_name, trait[key])
			return database
		except RuntimeError as error:
			print("Error with adding input (key {0}, trait {1}, trait_name {2}".format(key, trait, trait_name))
			open('../meta/error_log', 'a+').write(str(time.time())+ "\n\n"+ str(error))
			raise