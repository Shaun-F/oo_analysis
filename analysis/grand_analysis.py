"""
grand_analysis.py: generates grand spectraum for a number of significant figures

Created By: Erik Lentz
"""
import sys
sys.path.append("..")
import analysis.scan_analysis
import analysis.coaddition
import h5py
# seed result grand spectra

# add scan results to grand spectra

file = h5py.File(b"../data/raw/run1a_data.hdf5", "r+")

#pull grand spectra or create it. Then pull data group
if "grand_spectra_run1a" in file.keys():
	grand_spectra = file["grand_spectra_run1a"]
else:
	grand_spectra = file.create_dataset("grand_spectra_run1a", data = [], dtype=float, chunks = True, maxshape=None)

data_group = file["digitizer_log_run1a"]

# cycle over scans, analyzing each
for dataset in data_group.keys():
	data = dataset[...]












