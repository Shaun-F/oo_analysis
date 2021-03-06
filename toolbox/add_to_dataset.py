import numpy 
import numba; from numba import jit, njit


@jit
def addtodataset(dataset, array_or_string, position = None):
	"""
	First resize the dataset to accommodate more data then assign values to resized dataset
	position argument determines the index at which the array or string will be placed. Must be positive integer.
	"""
	ds = dataset[...] #copy dataset data into new array
	add = array_or_string

	if not (isinstance(add, str) and isinstance(add, bytes)):
		if position==None:
			ds = numpy.insert(ds, len(ds), add)
		else:
			ds = numpy.insert(ds, position, add)
	else:
		if position==None:
			ds = numpy.insert(ds, len(ds), add.encode()) #append to end of array
		else:
			ds = numpy.insert(ds, position, add.encode()) #insert into array at position
	#return modified data to dataset
	dataset.resize(ds.shape)
	dataset[...]=ds
	

def subtractfromdataset(dataset, array_or_string=None, position=None):
	ds = dataset[...]
	sub = array_or_string
	pos = position
	try:
		if array_or_string!=None:
			if not (isinstance(sub, str) or isinstance(sub, bytes)):
				for i in sub:
					ds = numpy.delete(ds, numpy.argwhere(ds==i))
			else:
				sub = sub.encode()
				ds = numpy.delete(ds, numpy.argwhere(ds==sub))
		elif position!=None:
			ds = numpy.delete(ds, [position])
		else:
			return "Error: Please specify an item to remove or an index"
		dataset.resize(ds.shape)
		dataset[...]=ds
	except IndexError:
		raise
		

def assign_newdata(dataset, data):
	"""
	Function reassigns the data array in a h5py dataset.
	NOTE: This overwrites all data in a dataset with the new data
	"""
	data = numpy.asarray(data)
	try:
		dataset.resize(data.shape)
		dataset[...] = data
	
	except (IndexError, Keyerror, MemoryError):
		raise
	
	
	
	
	
	