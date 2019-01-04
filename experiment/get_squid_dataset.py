
import sys
sys.path.append("..")

import h5py
import numpy as np
import dateutil.parser
from dateutil.parser import parse
def get_squid_dataset(timestamp):
	"""function gets the squid_temperature value closest to the input timestamp. Must be in format dd-mm-yyyy hh:mm:ss"""
	f = h5py.File(b"../data/raw/run1a_data.hdf5", "r")["squid_temperature_run1a"]
	keys = [x for x in f.keys()]
	dateobj = parse(timestamp, dateutil.parser.parserinfo(dayfirst=True)).strftime(("%d-%m-%Y %H:%M:%S"))
	
	year, month, day = dateobj[6:10], dateobj[3:5],dateobj[0:2]
	hour, minute, second = dateobj[11:13], dateobj[14:16], dateobj[17:19]
	
	###
	#This isnt a very robust retrieval of dataset, but I could not think of anything other way.
	###
	closest = []
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
							closest.append(x)
	if len(closest)==1:
		return closest[0]
	elif len(closest)!=0:
		minutes = []
		minute = int(minute)
		for time in closest:
			minutes.append(int(time[14:16])) #should have same indices as 'closest' list
		val, inx = min((abs(val), inx) for (inx, val) in enumerate([x-minute for x in minutes]))
		return f[closest[inx]]

		
			
