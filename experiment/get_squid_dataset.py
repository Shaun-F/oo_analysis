import h5py
import numpy as np
def get_squid_dataset(timestamp):
	"""function gets the squid_temperature value closest to the input timestamp. Must be in format dd-mm-yyyy hh:mm:ss"""
	f = h5py.File(u'../data/raw/run1a_data.hdf5', 'r')["squid_temperature_run1a"]
	keys = [x for x in f.keys()]
	year, month, day = str(timestamp[6:10]),str(timestamp[3:5]),str(timestamp[0:2])
	hour, minute, second = str(timestamp[11:13]),str(timestamp[14:16]),str(timestamp[17:19])
	#This isnt a very robust retrieval of dataset, but I could not think of anything other way.
	for x in keys:
		if (x[6:10]==year and x[3:5]==month and x[0:2]==day):
			if x[11:13]==hour:
				if x[14:16]==minute:
					return f[x]
				else:
					for i in np.arange(15):
						minute = int(minute) #turn into number
						if (x[14:16]==str(minute-i) or x[14:16]==str(minute+i)):
							#Find closest value
							return f[x]