# may now be defunct as I have tried to reproduce in core
# please make sure core method is working properly before deleting this code

import h5py
from get_squid_dataset import get_squid_dataset
from calc_sys_temp_offline import *

def determine_bad_scans(h5py_file=h5py.File('run1a_data', 'r+'), scan_number=0):
	"""Function takes a scan and decides if it should be cut from the grand spectrum.
	Current cut reasons: Timestamp, anomalous Q, Sys temp too high, no associated axion log, non-empty notes except for filled_in_after_sag
	"""
	cut = False
	cut_reason = ""
	#get relevant dataset parameters
	id_string = str(scan_number)
	dig_dataset = h5py_file["digitizer_log_run1a"][id_string]
	axion_dataset = h5py_file["axion_log_run1a"][id_string]
	timestamp = axion_dataset.attrs["timestamp"]
	bad_timestamps = h5py_file["bad_timestamps"]
	squid_dataset = get_squid_dataset(timestamp)
	Q = float(axion_dataset.attrs["Q"])
	squid_temp=float(squid_dataset[...])
	notes = axion_dataset.attrs["notes"]
	start = float(dig_dataset.attrs["start_frequency"])
	stop = float(dig_dataset.attrs["stop_frequency"])
	Tsys = calc_sys_temp(dig_dataset, axion_dataset, squid_dataset)["system temperature"]

	if timestamp in bad_timestamps:
		cut=True
		cut_reason = "Scan not suitable for analysis (Bad timestamp)"
	if Q>10**6:
		cut=True
		cut_reason = "Q too high"
	elif start<644.9 or stop<644.9 or start>685.1 or stop>685.1:
		cut=True
		cut_reason = "scan outside 645 to 685 MHz range"
	elif Tsys>5:
		cut=True
		cut_reason = "System temperature over 5k"
	elif notes!="nan" and notes!="filled_in_after_SAG":
		cut=True
		cut_reason = "Scan not suitable for analysis (See dataset notes)"
	else:
		#Sanity check
		cut = False
		cut_reason = ""

	axion_dataset.attrs["cut"] = cut
	axion_dataset.attrs["cut_reason"] = cut_reason
	dig_dataset.attrs["cut"] = cut
	dig_dataset.attrs["cut_reason"] = cut_reason


determine_bad_scans(scan_number=516305)
