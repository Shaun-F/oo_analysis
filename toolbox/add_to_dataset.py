import numpy 

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
	dataset.resize(ds.shape[0], axis=0)
	dataset[...]=ds
	
def subtractfromdataset(dataset, array_or_string=None, position=None):
	ds = dataset[...]
	sub = array_or_string
	pos = position
	try:
		if array_or_string!=None:
			if not (isinstance(sub, str) or isinstance(sub, bytes)):
				ds = numpy.delete(ds, numpy.argwhere(ds==sub))
			else:
				sub = sub.encode()
				ds = numpy.delete(ds, numpy.argwhere(ds==sub))
		elif position!=None:
			ds = numpy.delete(ds, [position])
		else:
			return "Error: Please specify an item to remove or an index"
		dataset.resize(ds.shape[0], axis=0)
		dataset[...]=ds
	except IndexError:
		raise
	
	
	
	
	
	