"""
Script plots the figures of merit for the analysis
"""
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os, sys
import pickle
import scipy.stats as stats
import copy
import datetime as dt
from datetime import timedelta as td
import dateutil
from dateutil.parser import parse
import calendar

"""
NOTE:
To Retrieve data from an AxesSubPlot object, use AxesSubPlot.get_lines()[0].get_ydata()
"""


class figures_class():
	
	def __init__(self, savedir='here', **kwargs):
		
		if savedir == 'here':
			savedir = os.getcwd() + "/oo_analysis/figures/"
		#set values common to all plots
		if 'object' not in kwargs.keys():
			raw_data_filename = os.getcwd() + "/oo_analysis/data/raw/run1a_data.hdf5"
			self.file = h5py.File(raw_data_filename.encode(), 'r')
			self.grand_spectra_group = self.file['grand_spectra_run1a']
			self.digitizer_group = self.file['digitizer_log_run1a']
			self.axion_frequencies_MHz = self.grand_spectra_group['axion_frequencies'][...]*10**(-6)
		else:
			self.object = kwargs['object']
			self.file = self.object.h5py_file
			self.grand_spectra_group = self.file['grand_spectra_run1a']
			self.digitizer_group = self.object.dig_dataset
			self.axion_frequencies_MHz = self.grand_spectra_group['axion_frequencies'][...]*10**(-6)
			
		#Plot parameters
		self.style = 'fast'
		plt.style.use(self.style)
		self.extension = ".pdf"
		self.savefile = savedir + '/'
		self.dpi = 1200
		self.label_font_size = 18
		self.legend_font_size = 11
		self.linewidth = 3.0
		self.alpha = 0.5
		self.xpad = 12
		self.show = False
		self.coupling_ylims = (0, 2)
		self.RF_interference_mask = (self.axion_frequencies_MHz<660.16)|(self.axion_frequencies_MHz>660.27)
		
		#force attributes from kwargs, overriding above
		for key,value in kwargs.items():
			setattr(self, key, value)
			
			
	def sensitivity_coupling(self, window_average = 0.1, **kwargs):
	
		data = self.grand_spectra_group['sensitivity_coupling'][...]
		mask = np.isfinite(data)
		master_mask = mask&self.RF_interference_mask
		reduced_data = data[master_mask]
		domain = self.axion_frequencies_MHz[master_mask]
		
		nbins = int((max(domain)-min(domain))/window_average)
		binner = stats.binned_statistic(domain, reduced_data, bins = nbins)
		reduced_data = binner.statistic
		domain = binner.bin_edges[:nbins]
		
		ax = plt.subplot(111)
		
		plt.title("Coupling sensitivity of run1A")
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel(r"Coupling Sensitivity ($g_{a \gamma \gamma}/g_{a \gamma \gamma, DFSZ}$)")
		
		plt.ylim(self.coupling_ylims[0], self.coupling_ylims[1])
		plt.locator_params(axis='x', nbins=8)
		plt.xticks(rotation=25)
		plt.ticklabel_format(useOffset=False)
		
		#plot
		plt.tight_layout()
		
		plt.plot(domain, reduced_data, ls='steps')
		plt.savefig(self.savefile + "Coupling_sensitivity" + self.extension, dpi = self.dpi) #Save image
		#Save plot object to file
		if 'params' in kwargs:
			params = kwargs['params']
			pickle.dump(ax, open(os.getcwd() + '/oo_analysis/figures/pickles/Coupling_sensitivity_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
		else:
			pickle.dump(ax, open(os.getcwd() + "/oo_analysis/figures/pickles/Coupling_sensitivity.pickle", 'wb'))
			
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		plt.clf()
	def sensitivity_power(self, window_average = 0.1, **kwargs):
		"""
		Plot the sensitivity to power"
		"""
		data = self.grand_spectra_group['sensitivity_power'][...]
		mask = np.isfinite(data)
		master_mask = mask&self.RF_interference_mask
		reduced_data = data[master_mask]
		domain = self.axion_frequencies_MHz[master_mask]
		
		
		nbins = int((max(domain)-min(domain))/window_average)
		binner = stats.binned_statistic(domain, reduced_data, bins = nbins)
		reduced_data = binner.statistic
		domain = binner.bin_edges[:nbins]
		
		ax = plt.subplot(111)
		plt.locator_params(axis='x', nbins=8)
		plt.xticks(rotation=25)
		plt.ticklabel_format(useOffset=False)

		plt.title("Power Sensitivity of Run1A")
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel(r"Power Sensitivity ($P/P_{DFSZ}$)")
		
		plt.ylim(self.coupling_ylims[0], self.coupling_ylims[1])

		plt.plot(domain, reduced_data, ls='steps')
		
		#plot
		plt.tight_layout()
		plt.savefig(self.savefile + "power_sensitivity" + self.extension, dpi = self.dpi)

		if 'params' in kwargs:
			params = kwargs['params']
			pickle.dump(ax, open(os.getcwd() + '/oo_analysis/figures/pickles/Power_sensitivity_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
		else:
			pickle.dump(ax, open(os.getcwd() + "/oo_analysis/figures/pickles/Power_sensitivity.pickle", 'wb'))
		
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		plt.clf()
	def SNR(self, **kwargs):
		"""
		Plot the SNR
		"""
		data = self.grand_spectra_group['SNR'][...]
		mask = np.isfinite(data)
		master_mask = mask&self.RF_interference_mask
		reduced_data = data[master_mask]
		domain = self.axion_frequencies_MHz[master_mask]
		
		ax = plt.subplot(111)
		plt.title("SNR of run1a")
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel("SNR")
		
		
		#plot
		plt.tight_layout()
		
		plt.plot(domain, reduced_data)
		plt.savefig(self.savefile + "SNR" + self.extension, dpi = self.dpi)
		
		if 'params' in kwargs:
			params = kwargs['params']
			pickle.dump(ax, open(os.getcwd() + '/oo_analysis/figures/pickles/SNR_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
		else:
			pickle.dump(ax, open(os.getcwd() + "/oo_analysis/figures/pickles/SNR.pickle", 'wb'))
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		plt.clf()
	def sensitivity_DM(self, window_average = 0.1, **kwargs):
		"""
		Plot the sensitivity to dark matter densities
		"""
		data = self.grand_spectra_group['sensitivity_power'][...]*0.45
		mask = np.isfinite(data)
		master_mask = mask&self.RF_interference_mask
		reduced_data = data[master_mask]
		domain = self.axion_frequencies_MHz[master_mask]
		
		nbins = int((max(domain)-min(domain))/window_average)
		binner = stats.binned_statistic(domain, reduced_data, bins = nbins)
		reduced_data = binner.statistic
		domain = binner.bin_edges[:nbins]
		
		
		ax = plt.subplot(111)
		
		plt.title("Dark Matter Density sensitivity of run1A")
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel("Dark Matter Density Sensitivity (DFSZ/GeV/cc)")
		
		plt.ylim(self.coupling_ylims[0], self.coupling_ylims[1])
		
		plt.tight_layout()
		plt.plot(domain, reduced_data, ls='steps')
		plt.savefig(self.savefile + "DM_sensitivity" + self.extension, dpi = self.dpi)
		
		if 'params' in kwargs:
			params = kwargs['params']
			pickle.dump(ax, open(os.getcwd() + '/oo_analysis/figures/pickles/DM_sensitivity_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
		else:
			pickle.dump(ax, open(os.getcwd() + "/oo_analysis/figures/pickles/DM_sensitivity.pickle", 'wb'))
		
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		plt.clf()	
	def fit_significance(self, hist=False, yscale='log', xscale='linear', **kwargs):
		"""
		Plot the axion fit significance
		"""
		if hist:
			self.data = self.grand_spectra_group['optimal_weight_sum'][...]/self.grand_spectra_group['model_excess_sqrd'][...]/self.grand_spectra_group['axion_fit_uncertainty'][...]
			mask = np.isfinite(self.data)
			self.master_mask = mask&self.RF_interference_mask
			reduced_data = self.data[self.master_mask]
			domain = self.axion_frequencies_MHz[self.master_mask]
			mean = np.mean(reduced_data)
			disp = np.std(reduced_data)
			
			ax = plt.subplot(111)
			plt.title("Axion Fit Significance Distribution")
			plt.xlabel("Normalized Power Excess")
			plt.ylabel("Count")
			
			plt.tight_layout()
			n, bins, patches = plt.hist(reduced_data, bins=10**3, align='mid', histtype='step', label='Fit Significance Distribution')
			plt.text(0, max(n)/15, r"$\mu$ = {0:0.2f}".format(mean) + "\n" + "$\sigma$ = {0:0.2f}".format(disp), horizontalalignment='center', verticalalignment='center')
			
			fit = stats.norm.pdf(bins, mean, disp)
			fit_plot = plt.plot(bins, fit*(max(n)/max(fit)), 'r--', lw=1, alpha=0.5, color='black',  label='Fitted Quadratic Monomial')
			
			plt.yscale(yscale)
			plt.xscale(xscale)
			plt.xlim(-7, 7)
			plt.ylim(10**(-1), 4*max(n))
			plt.legend()
			plt.savefig(self.savefile + "Fit_Significance_Distribution" + self.extension, dpi = self.dpi)
			
			if 'show' in list(kwargs.keys()) and kwargs['show']:
				plt.show()
			
			if 'params' in kwargs:
				params = kwargs['params']
				pickle.dump(ax, open(os.getcwd() + '/oo_analysis/figures/pickles/fit_significance_distribution_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
			else:
				pickle.dump(ax, open(os.getcwd() + "/oo_analysis/figures/pickles/fit_significance_distribution.pickle", 'wb'))
			plt.clf()
			
			
		else:
			self.data = self.grand_spectra_group['optimal_weight_sum'][...]/self.grand_spectra_group['model_excess_sqrd'][...]/self.grand_spectra_group['axion_fit_uncertainty'][...]
			mask = np.isfinite(self.data)
			mask1 = (self.data>0)
			self.master_mask = mask&self.RF_interference_mask&mask1
			reduced_data = self.data[self.master_mask]
			domain = self.axion_frequencies_MHz[self.master_mask]
			
			"""
			nbins = int(np.round(len(domain)/7)) #7 bins wide
			binner = stats.binned_statistic(domain, reduced_data, bins = nbins, statistic='max')
			reduced_data = binner.statistic
			domain = binner.bin_edges[:nbins]
			"""
			
			ax = plt.subplot(111)
			plt.title("Axion Fit Significance")
			plt.xlabel("Axion Frequencies (MHz)")
			plt.ylabel("Significance")
			
			plt.tight_layout()
			plt.plot(domain, reduced_data)
			plt.savefig(self.savefile + "Fit_Significance" + self.extension, dpi = self.dpi)
			
			if 'show' in list(kwargs.keys()) and kwargs['show']:
				plt.show()
			
			if 'params' in kwargs:
				params = kwargs['params']
				pickle.dump(ax, open(os.getcwd() + '/oo_analysis/figures/pickles/fit_significance_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
			else:
				pickle.dump(ax, open(os.getcwd() + "/oo_analysis/figures/pickles/fit_significance.pickle", 'wb'))
			plt.clf()
	def axion_fit_uncertainty(self, **kwargs):
		data = self.grand_spectra_group['axion_fit_uncertainty'][...]
		mask = np.isfinite(data)
		master_mask = mask&self.RF_interference_mask
		reduced_data = data[master_mask]
		domain = self.axion_frequencies_MHz[master_mask]
		
		
		ax = plt.subplot(111)
		plt.title("Axion Fit Uncertainty")
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel(r"Uncertainty")
		
		plt.tight_layout()
		plt.plot(domain, reduced_data)
		plt.savefig(self.savefile + "Axion_fit_Uncertainty" + self.extension, dpi = self.dpi)
		if 'show' in list(kwargs.keys()) and kwargs['show']:
			plt.show()
			
		if 'params' in kwargs:
			params = kwargs['params']
			pickle.dump(ax, open(os.getcwd() + '/oo_analysis/figures/pickles/fit_uncertainty_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
		else:
			pickle.dump(ax, open(os.getcwd() + "/oo_analysis/figures/pickles/fit_uncertainty.pickle", 'wb'))
			
		plt.clf()
	def time_v_freq(self, **kwargs):
		"""
		Plot the frequency scanned over time
		"""
		if 'timestamp' not in kwargs.keys():
			timestamp = [self.digitizer_group[key].attrs['timestamp'] for key in self.digitizer_group]
		else:
			timestamp = kwargs['timestamp']
		if 'mode_frequencies' not in kwargs.keys():
			mode_freqs = [float(copy.deepcopy(self.digitizer_group[key].attrs['mode_frequency'])) for key in self.digitizer_group]
		else:
			mode_freqs = kwargs['mode_frequencies']
		bad_scans_mask = [self.digitizer_group[key].attrs['cut'] for key in self.digitizer_group]
		good_scans_mask = np.invert(bad_scans_mask)
		
		times = list(map(lambda x: dt.datetime.strptime(x, "%d-%m-%Y %H:%M:%S"), times)) #set all timestamp strings to datetime objects
		
		fig, ax = plt.subplots()
		
		ax.scatter(times[bad_scans_mask], mode_freq[bad_scans_mask], color='red', s=1)
		ax.scatter(times[good_scans_mask], mode_freq[bad_scans_mask], color='blue', s=1)
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
		
		if 'show' in list(kwargs.keys()) and kwargs['show']:
			plt.show()
		
		plt.clf()		
	def scan_rate(self, **kwargs):
		"""
		Plot the number of scans per frequency   ### NOTE: Use logarithmic y-axis
		"""
		import dateutil
		from dateutil.parser import parse

		if 'object' not in kwargs.keys():
			self.bad_scans = [self.digitizer_group[key] for key in self.digitizer_group if 'cut' in self.digitizer_group[key].attrs and self.digitizer_group[key].attrs['cut']]
			self.timestamp = [self.digitizer_group[key].attrs['timestamp'] for key in self.digitizer_group if 'alog_timestamp' in self.digitizer_group[key].attrs and 'cut' in self.digitizer_group[key].attrs and not self.digitizer_group[key].attrs['cut']]
			self.mode_freqs = [float(copy.deepcopy(self.digitizer_group[key].attrs['mode_frequency'])) for key in self.digitizer_group if 'alog_timestamp' in self.digitizer_group[key].attrs and 'cut' in self.digitizer_group[key].attrs and not self.digitizer_group[key].attrs['cut']]
			self.Q = np.array([float(copy.deepcopy(self.digitizer_group[key].attrs['Q'])) for key in self.digitizer_group if 'alog_timestamp' in self.digitizer_group[key].attrs and 'cut' in self.digitizer_group[key].attrs and not self.digitizer_group[key].attrs['cut']])
			self.bad_timestamps = [parse(i.attrs['timestamp'], dayfirst=True) for i in self.bad_scans]
		else:
			self.bad_scans = list(self.object.bad_scans.values())
			self.bad_timestamps = [parse(i.attrs['timestamp'], dayfirst=True) for i in self.bad_scans]
			self.timestamp = list(self.object.timestamp.values())
			self.mode_freqs = list(self.object.mode_frequencies.values())
			self.Q = list(self.object.Q.values())
		
		self.timestamps = np.array([parse(i, dayfirst=True) for i in self.timestamp])
	
		timestamps_gm = [calendar.timegm(i.timetuple()) for i in self.timestamps]
		bad_timestamps_gm = [calendar.timegm(i.timetuple()) for i in self.bad_timestamps]
		
		Q_percent = (self.Q/70000*10)#Percentage relative to upper bound placed on analysis
		self.bad_yvalues = np.interp(bad_timestamps_gm, timestamps_gm, self.mode_freqs)
		self.bad_qvalues = np.interp(bad_timestamps_gm, timestamps_gm, self.Q)
		
		fig = plt.figure()
		plt.title("Mode Frequency versus Time of scan")
		plt.xlabel("Timestamp")
		plt.ylabel(r"Mode Frequency (MHz)")
		
		plt.locator_params(axis='x', nbins=8)
		plt.xticks(rotation=25)
		plt.fill_between(self.timestamps, self.mode_freqs-self.Q/70000*10, self.mode_freqs+self.Q/70000*10, alpha=0.4, color='green', label='Q Width relative to 7000')
		plt.plot(self.timestamps, self.mode_freqs, label='Mode Frequency', color='black')
		
		[plt.vlines(self.bad_timestamps[i], (self.bad_yvalues-self.bad_qvalues/7000)[i], (self.bad_yvalues + self.bad_qvalues/7000)[i], color=[1,0.3, 0.3, 0.3]) for i in range(len(self.bad_yvalues))]
		red_patch = mpatches.Patch(color='red', label='Bad Scans')
		
		
		plt.tight_layout()
		plt.legend()
		plt.savefig(self.savefile + "scan_rate" + self.extension, dpi=self.dpi)
		
		if 'show' in list(kwargs.keys()) and kwargs['show']:
			plt.show()
		
		if 'params' in kwargs:
			params = kwargs['params']
			pickle.dump(fig, open(os.getcwd() + '/oo_analysis/figures/pickles/mode_freq_v_scantime_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
		else:
			pickle.dump(fig, open(os.getcwd() + "/oo_analysis/figures/pickles/mode_freq_v_scantime.pickle", 'wb'))
		
		plt.clf()	
	def deltas(self, **kwargs):
		"""
		Plot the deltas as histogram. Measures gaussianity of data. ### stack all deltas on top of one another and measure dispersion that way + distribution.
		"""
		deltas = self.grand_spectra_group['power_deviation'][...]
		sigma = np.std(deltas)
		nsigma = 3
		mask = (deltas>(-nsigma*sigma))|(deltas<(nsigma*sigma))
		deltas = deltas[mask]/sigma
		ax = plt.subplot(111)
		n, bins,patch = plt.hist(deltas, bins=10**4, align='mid', histtype='step')
		mean = np.abs(np.mean(deltas))
		disp = np.std(deltas)
		plt.yscale('log')
		plt.xscale('linear')
		plt.ylabel("Count")
		plt.xlabel("Normalized Power Excess")
		plt.text(-.2, max(n)/15, r"$\mu$ = {0:0.2f}".format(mean) + "\n" + "$\sigma$ = {1:0.2f}".format(mean, disp))
		plt.xlim(-nsigma, nsigma)
		plt.savefig(self.savefile + 'Run1a_deltas' + self.extension, dpi = self.dpi)
		
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		
		if 'params' in kwargs:
			params = kwargs['params']
			pickle.dump(ax, open(os.getcwd() + '/oo_analysis/figures/pickles/deltas_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
		else:
			pickle.dump(ax, open(os.getcwd() + "/oo_analysis/figures/pickles/deltas.pickle", 'wb'))
		
		
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
			timestamps = np.array([parse(i, dayfirst=True) for i in kwargs['timestamp']])
		else:
			for key in self.digitizer_group:
				attrs = self.digitizer_group[key].attrs
				if 'alog_timestamp' in attrs and not attrs['cut']:
					temperatures.append(attrs['squid_temperature'])
					timestamps.append(parse(attrs['timestamp'], dayfirst=True))
				
		fig, ax = plt.subplots()
		plt.tight_layout()
		plt.plot(timestamps, temperatures, lw=1)
		
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
	
		