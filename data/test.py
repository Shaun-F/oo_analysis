
test = "true"

def test():
	import h5py
	import os
	print(os.getcwd())
	print(h5py.File(u"raw/run1a_data.hdf5", "r+"))
	
test()