"""
Script plots the figures of merit for the analysis
"""
import h5py
import matplotlib.pyplot as plt
plt.style.use('seaborn-darkgrid')
import numpy as np
import os, sys



class figures_class():
	
	def __init__(self, savedir='here', **kwargs):
		
		if savedir == 'here':
			savedir = os.path.abspath(os.path.dirname(sys.argv[0])) + "/oo_analysis/figures/"
		#set values common to all plots
		raw_data_filename = os.getcwd() + "/oo_analysis/data/raw/run1a_data.hdf5"
		self.file = h5py.File(raw_data_filename.encode(), 'r')
		self.grand_spectra_group = self.file['grand_spectra_run1a']
		self.digitizer_group = self.file['digitizer_log_run1a']
		self.axion_frequencies_MHz = self.grand_spectra_group['axion_frequencies'][...]*10**(-6)
		
		#Plot parameters
		self.extension = ".pdf"
		self.savefile = savedir + '\\'
		self.dpi = 600
		self.label_font_size = 18
		self.legend_font_size = 11
		self.linewidth = 3.0
		self.alpha = 0.5
		self.xpad = 12
		self.show = False
		self.RF_interference_mask = (self.axion_frequencies_MHz<660.16)|(self.axion_frequencies_MHz>660.27)
		
		#force attributes from kwargs, overriding above
		for key,value in kwargs.items():
			setattr(self, key, value)
			
			
	def sensitivity_coupling(self, **kwargs):
		data = self.grand_spectra_group['sensitivity_coupling'][...]
		mask = np.isfinite(data)
		master_mask = mask&self.RF_interference_mask
		reduced_data = data[master_mask]
		domain = self.axion_frequencies_MHz[master_mask]
		plt.title("Coupling sensitivity of run1a")
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel(r"Coupling Sensitivity (1/DFSZ)")
		plt.locator_params(axis='x', nbins=8)
		plt.xticks(rotation=25)
		plt.ticklabel_format(useOffset=False)
		
		#plot
		plt.tight_layout()
		
		plt.plot(domain, reduced_data)
		plt.savefig(self.savefile + "Coupling_sensitivity" + self.extension, dpi = self.dpi)
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		plt.clf()
	def sensitivity_power(self, **kwargs):
		"""
		Plot the sensitivity to power"
		"""
		data = self.grand_spectra_group['sensitivity_power'][...]
		mask = np.isfinite(data)
		master_mask = mask&self.RF_interference_mask
		reduced_data = data[master_mask]
		domain = self.axion_frequencies_MHz[master_mask]
		plt.locator_params(axis='x', nbins=8)
		plt.xticks(rotation=25)
		plt.ticklabel_format(useOffset=False)

		plt.title("power sensitivity of run1a")
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel(r"Power Sensitivity")
		
		
		#plot
		plt.tight_layout()
		
		plt.plot(domain, reduced_data)
		plt.savefig(self.savefile + "power_sensitivity" + self.extension, dpi = self.dpi)
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		plt.clf()
	def SNR(self, **kwargs):
		"""
		Plot the SNR
		"""
		data = self.grand_spectra_group['SNR'][...]
		mask = np.isfinite(data)
		reduced_data = data[mask]
		domain = self.axion_frequencies_MHz[mask]
		plt.title("SNR of run1a")
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel("SNR")
		
		
		#plot
		plt.tight_layout()
		
		plt.plot(domain, reduced_data)
		plt.savefig(self.savefile + "SNR" + self.extension, dpi = self.dpi)
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		plt.clf()
	def sensitivity_DM(self, **kwargs):
		"""
		Plot the sensitivity to dark matter densities
		"""
		data = self.grand_spectra_group['sensitivity_power'][...]
		mask = np.isfinite(data)
		reduced_data = data[mask]
		domain = self.axion_frequencies_MHz[mask]
		plt.title("Sensitivity to Dark Matter density at DSFZ level")
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel("Dark Matter Density Sensitivity")
		
		plt.tight_layout()
		plt.plot(domain, reduced_data)
		plt.savefig(self.savefile + "DM_sensitivity" + self.extension, dpi = self.dpi)
		
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		plt.clf()	
	
	def axion_fit_uncertainty(self, **kwargs):
		data = self.grand_spectra_group['axion_fit_uncertainty'][...]
		mask = np.isfinite(data)
		reduced_data = data[mask]
		domain = self.axion_frequencies_MHz[mask]
		plt.title("Axion Fit Uncertainty")
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel(r"Uncertainty")
		
		plt.tight_layout()
		plt.plot(domain, reduced_data)
		plt.savefig(self.savefile + "Axion_fit_Uncertainty" + self.extension, dpi = self.dpi)
		if 'show' in list(kwargs.keys()) and kwargs['show']:
			plt.show()
		plt.clf()
		
		
	def time_v_freq(self, **kwargs):
		"""
		Plot the frequency scanned over time
		"""
		import datetime as dt
		mode_freq = []
		times = []
		for i in self.digitizer_group:
			if 'alog_timestamp' in self.digitizer_group[i].attrs and not self.digitizer_group[i].attrs['cut']: #restrict to scans that have an associated axion scan log
					mode_freq.append(self.digitizer_group[i].attrs['mode_frequency'])
					times.append(self.digitizer_group[i].attrs['timestamp'])

		times = list(map(lambda x: dt.datetime.strptime(x, "%d-%m-%Y %H:%M:%S"), times)) #set all timestamp strings to datetime objects
		
		fig, ax = plt.subplots()
		
		ax.scatter(times, mode_freq, s=2)
		ax.set_xlabel("Date")
		ax.set_ylabel("Mode Frequency (MHz)")
		ax.set_title("Mode Frequency versus date")
		
		#set ticks
		plt.locator_params(axis='x', nbins=8)
		plt.xticks(rotation=25)
		plt.tight_layout()
		plt.savefig(self.savefile + "time_v_freq" + self.extension, dpi = self.dpi)
		
		if 'show' in list(kwargs.keys()) and kwargs['show']:
			plt.show()
		plt.clf()	
	def candidates(self, **kwargs):
		"""
		Plot the placement of candidates over frequency
		"""
		#candidates as pulled from elog
		nib1candfreq = [652.37,652.43,652.45,652.505,652.545,655.74,655.76,655.965,656.00,656.04,656.08,658.54,659.48,659.565]
		nib1cand = [[freq,3.] for freq in nib1candfreq]
		nib2cand = [[ 660.19, 4.7780322225437], [ 660.2, 5.4971375979476], [ 660.24, 3.6087460437092], [ 660.245, 4.2776040976767], [ 660.295, 3.6891690650153], [ 660.315, 3.5036515112248], [ 660.325, 4.6918355682087], [ 660.82, 3.4991249763781], [ 660.915, 3.0200994610868], [ 660.925, 3.2868675329851], [ 660.96, 3.1940864636256], [ 660.985, 3.009667492588], [ 663.835, 6.54937474606], [ 663.845, 14.583801581512], [ 663.855, 4.6943727578531], [ 666.395, 3.0232789874931], [ 669.03, 7.7599761841003], [ 669.04, 9.4642824474134]]
		nib3cand = [[ 673.66, 3.1291448682523], [ 673.69, 5.5142109453939], [ 674.135, 3.3500074464926], [ 674.19, 3.8189878005964], [ 674.425, 4.5005720049264], [ 674.58, 4.0019504106299], [ 674.59, 3.9921224210076], [ 674.615, 3.1390728172432], [ 674.785, 5.1759340862907], [ 674.855, 3.6536106696123], [ 676.635, 3.2671759985871], [ 676.66, 5.736775837645], [ 676.67, 4.7424480481658], [ 676.68, 3.8134711888459], [ 676.69, 3.4157645909233], [ 676.715, 4.2304356137148], [ 676.725, 3.8698543797455]]
		nib4cand = []
		candlist = nib1cand+nib2cand+nib3cand+nib4cand
		candfreqs = [cand[0] for cand in candlist]
		
		# using contour, neet to creat backdrop
		freqs = np.linspace(645.,680.,(680-645)*200,endpoint=False)
		y = np.array([0,2])
		X, Y = np.meshgrid(freqs,y)
		cand = np.array([[2.,2.] if freq in candfreqs else [0.,0.] for freq in freqs]).T
			
			
		plt.clf()
		fig = plt.contourf(X,Y,cand)
		
		plt.yscale('linear')
		plt.xscale('linear')
		plt.xlabel(r'Frequency [MHz]', fontsize = self.label_font_size)
		
		plt.tick_params(axis='both', which='major', labelsize = self.label_font_size)
		plt.tick_params(axis='both', which='minor', labelsize = self.label_font_size)
		plt.axes().set_aspect('equal')
		plt.yticks([])
		
		#plt.legend(loc=1, prop={'size':self.legend_font_size})
		plt.tight_layout()
		plt.tick_params(axis='x', pad = self.xpad)

		plt.savefig(self.savefile + "Run1A_candidates" + self.extension, dpi = self.dpi)
		
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		
		plt.clf()		
	def scan_rate(self, **kwargs):
		"""
		Plot the number of scans per frequency   ### NOTE: Use logarithmic y-axis
		"""
		import matplotlib.image as img
		image = img.imread('scan_rate.png')
		plt.imshow(image)
		plt.axis('off')
		if 'show' in list(kwargs.keys()) and kwargs['show']:
			plt.show()
	
		plt.clf()	
	def deltas(self, **kwargs):
		"""
		Plot the deltas as histogram. Measures gaussianity of data. ### stack all deltas on top of one another and measure dispersion that way + distribution.
		"""
		deltas = self.grand_spectra_group['power_deviation'][...]
		plt.hist(deltas, bins=400)
		
		plt.savefig(self.savefile + 'Run1a_deltas' + self.extension, dpi = self.dpi)
		
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		
		plt.clf()	
	def injected_signals(self, **kwargs):
		"""
		Plot the injected signals. Confirmation of sensitivity. Gray might have this code......
		"""
	
		plt.clf()	
	def axion_model_compare(self, **kwargs):
		"""
			
		"""
	
		plt.clf()	
	def temp_v_freq(self, **kwargs):
		"""
		Plot the temperature over frequency
		
		
		
		Make sort of like the width-by-Q width plot. Since we have multiple passes, this will be multiple valued
		"""
		#set data
		temperatures = []
		frequencies = []
		keys = self.digitizer_group.keys()
		if 'Tsys' and 'mode_frequency' in list(kwargs.keys()):
			temperatures = kwargs['Tsys']
			frequencies = kwargs['mode_frequencies']
		else:
			for key in keys:
				attrs = self.digitizer_group[key].attrs
				if 'alog_timestamp' in attrs and not attrs['cut']:
					temperatures.append(attrs['squid_temperature'])
					frequencies.append(attrs['mode_frequency'])
		frequencies = np.asarray(frequencies)
		temperatures = np.asarray(temperatures)
		mask = (frequencies>0)
		frequencies = frequencies[mask]
		temperatures = temperatures[mask]
		plt.title("Temperature of system vs mode frequency")
		plt.xlabel("Frequency")
		plt.ylabel("Temperature (K)")
		plt.tight_layout()
		plt.scatter(frequencies, temperatures, s=2)
		plt.savefig(self.savefile + "Temperature_v_Frequency" + self.extension, dpi = self.dpi)
		
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		
		plt.clf()	
	def temp_v_time(self, **kwargs):
		"""
		Plot the temperature over time
		Format of digitizer log timestamps is %d-%m-%Y %H:%M:%S according to python datetime directives. 
		"""
		import datetime as dt
		
		
		#set data
		temperatures = []
		timestamps = []
		if 'Tsys' and 'timestamp' in list(kwargs.keys()):
			temperatures = kwargs['Tsys']
			timestamps = [dt.datetime.strptime(i, "%d-%m-%Y %H:%M:%S") for i in kwargs['timestamp']]
		else:
			for key in self.digitizer_group:
				attrs = self.digitizer_group[key].attrs
				if 'alog_timestamp' in attrs and not attrs['cut']:
					temperatures.append(attrs['squid_temperature'])
					timestamps.append(dt.datetime.strptime(attrs['timestamp'], "%d-%m-%Y %H:%M:%S"))
				
		fig, ax = plt.subplots()
		plt.tight_layout()
		ax.scatter(timestamps, temperatures, s=2)
		
		#set labels
		ax.set_title('Squid temperature vs data')
		ax.set_xlabel("Date")
		ax.set_ylabel("Temperature (K)")
		
		
		#set ticks
		plt.locator_params(axis='x', nbins=8)
		plt.xticks(rotation=25)
		
		plt.savefig(self.savefile + "temp_v_time" + self.extension, dpi = self.dpi)
		
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		plt.clf()	
	
		