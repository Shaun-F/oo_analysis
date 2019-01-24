"""
__init__.py: main init file for data management

Created by: Erik Lentz
Creation Date: 10/26/18
"""
import h5py
import pickle #used for serializing and de-serializing data.

def write_out(dataset,path):
	"""
	Function takes in a dataset, a file path and writes both the data and attributes to the file via serialization
	"""
	#data = h5py.File(path, 'w')
	# use write functioins to put grand spectra to file
	data = dataset[...]
	attributes = {key: dataset.attrs[key] for key in dataset.attrs}

	output_file = open(path, 'wb')
	pickle.dump(data, output_file)
	pickle.dump(attributes, output_file)

	output_file.close()
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