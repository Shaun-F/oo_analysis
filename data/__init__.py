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
	data = h5py.File(b"../data/raw/run1a_data.hdf5", "r+")
	
	if "bad_timestamps_run1a" not in data.keys():
		data.create_dataset(name="bad_timestamps_run1a", dtype="S10", data = [b"initval"], maxshape=(None,))
		
	dig_dataset = {}
	axion_dataset = {}
	squid_dataset = {}
	for i in range(start,stop):
		try:
			inx = str(i)
			dig_dataset[inx] = data["digitizer_log_run1a"][inx]
		except KeyError:
			continue #No digitizer log found. check next scan number
		try:
			inx = str(i)
			axion_dataset[inx] = data["axion_log_run1a"][inx]
		except KeyError:
			print(Exception("Error: No associated axion log for scan {0:s}".format(str(i))))
			continue
		try:
			inx=str(i)
			squid_dataset[inx] = get_squid_dataset(axion_dataset[inx].attrs["timestamp"])
		except KeyError:
			raise Exception("Error: No squid dataset found for scan {0:s}".format(str(i)))
	return dig_dataset, axion_dataset, squid_dataset, data

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
