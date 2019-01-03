"""
__init__.py: main init file for experiment management

Created by: Erik Lentz
Creation Date: 10/26/18
"""



import calc_sys_temp_offline
import get_squid_dataset

# get squid dataset
squid_dataset = get_squid_dataset.get_squid_dataset(timestamp) #need to make plural
squid_temp_dataset = float(squid_dataset[...]) # need to make plural
# calc sys temperature
def sys_temp(dig_data,axion_data): # need to make consistant with plural sets
    return calc_sys_temp_offline.calc_sys_temp(h5py_digitizer_dataset, h5py_axion_dataset, squid_temp_dataset)
