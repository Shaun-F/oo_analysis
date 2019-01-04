"""
__init__.py: main init file for experiment management

Created by: Erik Lentz
Creation Date: 10/26/18
"""


import sys
sys.path.append("..")
import experiment.calc_sys_temp_offline
import experiment.get_squid_dataset

# get squid dataset
def get_datasets(timestamps):
	if isinstance(timestamps, str):
		squid_dataset = get_squid_dataset.get_squid_dataset(timestamps)
		squid_temps = float(squid_dataset[...])
		
	elif isinstance(timestamps, list):
		squid_dataset = []
		squid_temps = []
		for timestamp in timestamps:
			dset = get_squid_dataset.get_squid_dataset(timestamp)
			squid_dataset.append(dset)
			squid_temps.append(float(dset[...]))
			
	elif isinstance(timestamps, dict):
		squid_dataset = {}
		squid_temps = {}
		for key, timestamp in timestamps:
			dset = get_squid_dataset.get_squid_dataset(timestamp)
			squid_dataset[key] = dset
			squid_temps[key] = float(dset[...])	
	else:
		return "Error: Timestamp container not valid. Must be type(str), type(list), or type(dict)"
	return squid_dataset, squid_temps
			
# calc sys temperature
def sys_temp(dig_data,axion_data): # need to make consistant with plural sets
    return calc_sys_temp_offline.calc_sys_temp(h5py_digitizer_dataset, h5py_axion_dataset, squid_temp_dataset)
