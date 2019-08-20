"""
__init__.py: main init file for data management

Created by: Erik Lentz
Creation Date: 10/26/18
"""
import h5py
import pickle #used for serializing and de-serializing data.

def write_out(datagroup,path):
	"""
	Function takes in a dataset, a file path and writes the data to the file via serialization
	"""
	#data = h5py.File(path, 'w')
	# use write functioins to put grand spectra to file
	
	datasets = {key: datagroup[key][...] for key in datagroup if key!="deltas"} #form dictionary of data to save to disk.
	output_file = open(path, 'wb') #file object for pickle to save to
	pickle.dump(datasets, output_file) #save dictionary to file object
	output_file.close() #close out file object
	
	return None
	
def read_file(file):
	"""
	De-serialize a file and read
	"""
	input = open(file, 'rb')

	output = []

	while True:
		try:
			output.append(pickle.load(input))
		except EOFError:
			break

	return output