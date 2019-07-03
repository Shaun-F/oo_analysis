import matplotlib.pyplot as plt
import numpy
import h5py
import warnings
warnings.filterwarnings("ignore")
def plotter(scan_number_or_array, savedir = None, show=False, **kwargs):
	#retrieve data
	dpi=1000
	plt.clf()
	whatplot = scan_number_or_array
	if isinstance(whatplot, int):
		scan_number = whatplot
		try:
			f = h5py.File(b"../data/raw/run1a_data.hdf5", "r")
			data_file = f['digitizer_log_run1a'][str(scan_number)]

			fstart = float(data_file.attrs['start_frequency'])
			fstop = float(data_file.attrs['stop_frequency'])
			fres = float(data_file.attrs['frequency_resolution'])
			
			data = data_file[...]
			data = data[numpy.isfinite(data)]
			
			
			domain = numpy.asarray([ fstart + fres*i for i in range(len(data))])
			
			if min(data)<0:
				m = min(data)*1.5
			elif min(data)>0:
				m=min(data)*0.75
			plt.ylim(m, max(data)*1.5)
			plt.title("Scan number {0}".format(scan_number))
			plt.xlabel("frequency (MHz)")
			plt.ylabel("power (normed)")
			for arg,val in kwargs.items():
				if arg=='title':
					plt.title(val)
					del kwargs['title']
				if arg=='xlabel':
					plt.xlabel(val)
					del kwargs['xlabel']
				if arg=='ylabel':
					plt.ylabel(val)
					del kwargs['ylabel']
				if arg=='ylimits':
					plt.ylim(val[0], val[1])
					del kwargs['ylimits']
				if arg=='xlimits':
					plt.xlim(val[0], val[1])
					del kwargs['xlimits']
				if arg=='xticks':
					plt.xticks(val)
					del kwargs['xticks']
				if arg=='yticks':
					plt.yticks(val)
					del kwargs['yticks']
			
			plt.locator_params(axis='x', nbins=12)
			plt.xticks(rotation=30)
			plt.tight_layout()
			plt.ticklabel_format(useOffset=False)
			plt.plot(domain, data, **kwargs)
			if savedir!=None:
				plt.savefig(savedir, dpi=dpi)
			if show:
				plt.show()
			f.close()
		except (AttributeError, KeyError):
			raise
	
	elif isinstance(whatplot, numpy.ndarray):
			if min(whatplot)<0:
				m = min(whatplot)
			elif min(whatplot)>=0:
				m=min(whatplot)
			ydata = whatplot[numpy.isfinite(whatplot.real)].real
			xdata = numpy.arange(len(ydata))
			plt.ylim(m, max(ydata)*(1+1/24))
			
			for arg,val in kwargs.items():
				if arg=='title':
					plt.title(val)
				if arg=='xlabel':
					plt.xlabel(val)
				if arg=='ylabel':
					plt.ylabel(val)
				if arg=='ylimits':
					plt.ylim(val[0], val[1])
				if arg=='xlimits':
					plt.xlim(val[0], val[1])
				if arg=='xticks':
					plt.xticks(val)
				if arg=='yticks':
					plt.yticks(val)
				if arg=='data':
					xdata = kwargs['data'][0]
					ydata = kwargs['data'][1]
			for string in ['title', 'xlabel', 'ylabel', 'xlimits', 'ylimits', 'xticks', 'yticks', 'data']:
				if string in kwargs.keys():
					del kwargs[string]
			
			plt.plot(xdata, ydata, **kwargs)
			
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
			for arg,val in kwargs.items():
				if arg=='title':
					plt.title(val)
				if arg=='xlabel':
					plt.xlabel(val)
				if arg=='ylabel':
					plt.ylabel(val)
				if arg=='ylimits':
					plt.ylim(val[0], val[1])
				if arg=='xlimits':
					plt.xlim(val[0], val[1])
				if arg=='xticks':
					plt.xticks(val)
				if arg=='yticks':
					plt.yticks(val)
				if arg=='data':
					xdata = kwargs['data'][0]
					ydata = kwargs['data'][1]
					
			plt.plot(numpy.arange(len(whatplot)), list(map(abs, whatplot)))

			if savedir!=None:
				plt.savefig(savedir, dpi=dpi)
			else:
				plt.show()
	
	plt.clf()