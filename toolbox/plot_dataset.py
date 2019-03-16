import matplotlib.pyplot as plt
import numpy
import h5py
def plotter(scan_number_or_array, savedir = None):
	#retrieve data
	dpi=1000
	
	whatplot = scan_number_or_array
	if isinstance(whatplot, int):
		scan_number = whatplot
		try:
			f = h5py.File(b"../data/raw/run1a_data.hdf5", "r")
			data_file = f['digitizer_log_run1a'][str(scan_number)]

			fstart = data_file.attrs['start_frequency']
			fstop = data_file.attrs['stop_frequency']
			fres = data_file.attrs['frequency_resolution']
			
			data = data_file[...]
			domain = numpy.asarray([ fstart + fres*i for i in range(len(data))])
			
			if min(data)<0:
				m = min(data)*1.5
			elif min(data)>0:
				m=min(data)*0.75
			plt.ylim(m, max(data)*1.5)
			plt.title("Scan number {0}".format(scan_number))
			plt.xlabel("frequency (hz)")
			plt.ylabel("power (normed)")
			
			
			plt.plot(domain, data)
			if savedir!=None:
				plt.savefig(savedir, dpi=dpi)
			else:
				plt.show()
			f.close()
		except (AttributeError, KeyError):
			raise
	elif isinstance(whatplot, numpy.ndarray):
			if min(whatplot)<0:
				m = min(whatplot)
			elif min(whatplot)>=0:
				m=min(whatplot)
			plt.ylim(m, max(whatplot))
			plt.plot(numpy.arange(len(whatplot)), whatplot.real)
			
			if savedir!=None:
				plt.savefig(savedir, dpi=dpi)
			else:
				plt.show()
	elif isinstance(whatplot, list):
			if min(whatplot)<0:
				m = min(whatplot)*1.5
			elif min(whatplot)>0:
				m=min(whatplot)*0.75
			plt.ylim(m, max(whatplot)*1.5)
			plt.plot(numpy.arange(len(whatplot)), list(map(abs, whatplot)))

			if savedir!=None:
				plt.savefig(savedir, dpi=dpi)
			else:
				plt.show()
	
	plt.clf()