"""
Script plots the figures of merit for the analysis

Created by: Shaun Fell
"""
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import StrMethodFormatter
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
To Retrieve data from an AxesSubPlot object, use AxesSubPlot.get_lines()[0].get_ydata() or AxesSubPlot.get_lines()[0].get_xdata()
"""


class figures_class():
	
	def __init__(self, savedir='here', import_hdf5 = True, **kwargs):
		
		if savedir == 'here':
			savedir = os.getcwd() + "/oo_analysis/figures/"
		
		#File path to paper
		paper_savedir = "D:/Users/shaun/Documents/ADMX/Papers/Own Papers/Orbital_Considerations_and_other_Advances_in_searching_for_Invisible_Axion_Dark_Matter_in_ADMX_Run1A/"
		self.save_to_paper = False
		if 'save_to_paper' in list(kwargs.keys()) and kwargs['save_to_paper']:
			if os.path.exists(paper_savedir):
				self.save_to_paper = kwargs['save_to_paper']
				self.paper_dir = paper_savedir
				
		
		
		#set values common to all plots
		if import_hdf5:
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
			
			self.RF_interference_mask = (self.axion_frequencies_MHz<660.16)|(self.axion_frequencies_MHz>660.27) #This frequency range contains interference from RF signals

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
		
		#force attributes from kwargs, overriding above
		for key,value in kwargs.items():
			setattr(self, key, value)
			
			
	def sensitivity_coupling(self, window_average = 1.0, **kwargs):
	
		data = self.grand_spectra_group['sensitivity_coupling'][...]
		mask = np.isfinite(data) #Mask out infinite values
		master_mask = mask&self.RF_interference_mask
		reduced_data = data[master_mask]
		domain = self.axion_frequencies_MHz[master_mask]
		
		#Bins coupling with bin size determined by window_average
		nbins = int(np.ceil((max(domain)-min(domain))/window_average))
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
		if self.save_to_paper:
			plt.savefig(self.paper_dir +  "Coupling_sensitivity" + self.extension, dpi = self.dpi) #Save image
			
		#Save plot object to file
		if 'params' in kwargs:
			params = kwargs['params']
			pickle.dump(ax, open(os.getcwd() + '/oo_analysis/figures/pickles/Coupling_sensitivity_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
		else:
			pickle.dump(ax, open(os.getcwd() + "/oo_analysis/figures/pickles/Coupling_sensitivity.pickle", 'wb'))
			
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		plt.clf()
	def sensitivity_power(self, window_average = 1.0, **kwargs):
		"""
		Plot the sensitivity to power"
		"""
		data = self.grand_spectra_group['sensitivity_power'][...]
		mask = np.isfinite(data)
		master_mask = mask&self.RF_interference_mask
		reduced_data = data[master_mask]
		domain = self.axion_frequencies_MHz[master_mask]
		
		#Bins coupling with bin size determined by window_average
		nbins = nbins = int(np.ceil((max(domain)-min(domain))/window_average))
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
		if self.save_to_paper:
				plt.savefig(self.paper_dir +  "power_sensitivity" + self.extension, dpi = self.dpi) #Save image
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
	def sensitivity_DM(self, window_average = 1.0, **kwargs):
		"""
		Plot the sensitivity to dark matter densities
		"""
		data = self.grand_spectra_group['sensitivity_power'][...]*0.45
		mask = np.isfinite(data)
		master_mask = mask&self.RF_interference_mask
		reduced_data = data[master_mask]
		domain = self.axion_frequencies_MHz[master_mask]
		
		#Bins coupling with bin size determined by window_average
		nbins = int(np.ceil((max(domain)-min(domain))/window_average))
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
		if self.save_to_paper:
			plt.savefig(self.paper_dir +  "DM_sensitivity" + self.extension, dpi = self.dpi) #Save image
		
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
			#For some reason, the axion fit significance doesnt coadd well. I never found the time to track it down, so this is my work around
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
			if self.save_to_paper:
				plt.savefig(self.paper_dir +  "Fit_Significance_Distribution" + self.extension, dpi = self.dpi) #Save image
			if 'show' in list(kwargs.keys()) and kwargs['show']:
				plt.show()
			
			if 'params' in kwargs:
				params = kwargs['params']
				pickle.dump(ax, open(os.getcwd() + '/oo_analysis/figures/pickles/fit_significance_distribution_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
			else:
				pickle.dump(ax, open(os.getcwd() + "/oo_analysis/figures/pickles/fit_significance_distribution.pickle", 'wb'))
			plt.clf()
			
			
		else:
			#For some reason, the axion fit significance doesnt coadd well. I never found the time to track it down, so this is my work around
			self.data = self.grand_spectra_group['optimal_weight_sum'][...]/self.grand_spectra_group['model_excess_sqrd'][...]/self.grand_spectra_group['axion_fit_uncertainty'][...]
			mask = np.isfinite(self.data)
			mask1 = (self.data>0)
			self.master_mask = mask&self.RF_interference_mask&mask1
			reduced_data = self.data[self.master_mask]
			domain = self.axion_frequencies_MHz[self.master_mask]
			
			if 'synthetic_injections' in kwargs.keys(): #Mask out synthetically injected axion signals
				synth_mask = [True]*len(reduced_data)
				for freq in kwargs['synthetic_injections']:
					synth_mask = synth_mask&(((domain*10**6<(float(freq)-1000)))|(domain*10**6>(float(freq)+1000)))
			
				reduced_data = reduced_data[synth_mask]
				domain = domain[synth_mask]
			
			
			#Bin 3-sigma excursions together that are close by
			three_sig_exc = reduced_data[reduced_data>3]
			three_sig_exc_domain = domain[reduced_data>3]
			
			
			res = domain[1]-domain[0]
			bins = np.arange(min(domain), max(domain), 7*0.001)
			indices = np.digitize(three_sig_exc_domain, bins)
			self.candidates_freq = [three_sig_exc_domain[indices==i].min() for i in np.unique(indices)]
			self.candidates = [three_sig_exc[indices==i].max() for i in np.unique(indices)]
						
			iter = 0
			self.binned_exc_dom = []
			self.binned_exc = []
			while iter<len(three_sig_exc_domain):
				seed = three_sig_exc_domain[iter]
				group = np.where(((-7*0.0000954 + seed) < three_sig_exc_domain)&(three_sig_exc_domain<(7*0.0000954 + seed)))
				middle_freq = three_sig_exc_domain[group][int(np.floor(len(three_sig_exc_domain[group])/2))]
				self.binned_exc_dom.append(middle_freq)
				
				self.binned_exc.append(max(three_sig_exc[group]))

				iter += len(group)


			ax = plt.subplot(111)
			plt.title("Axion Fit Significance")
			plt.xlabel("Axion Frequencies (MHz)")
			plt.ylabel("Significance")
			
			plt.tight_layout()
			plt.plot(domain, reduced_data, label='Fit Significance')
			plt.scatter(self.binned_exc_dom, self.binned_exc, color='black', s=5, label='Candidates ({0})'.format(len(self.binned_exc)))
			plt.axhline(3, label=r"3-$\sigma$", color='black', linestyle='--', alpha=0.7)
			plt.legend(loc='upper right')
			plt.savefig(self.savefile + "Fit_Significance" + self.extension, dpi = self.dpi)
			if self.save_to_paper:
				plt.savefig(self.paper_dir +  "Fit_Significance" + self.extension, dpi = self.dpi) #Save image
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
		
		#Mask out bad scans
		bad_scans_mask = [self.digitizer_group[key].attrs['cut'] for key in self.digitizer_group]
		good_scans_mask = np.invert(bad_scans_mask)
		
		#Generate list of datetime objects
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
		self.candlist = nib1cand+nib2cand+nib3cand+nib4cand #49 candidates
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
		Plot the mode frequency versus timestamp. Width of line determined by Q factor
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
	def deltas(self, hist=False, **kwargs):
		"""
		Plot the deltas as histogram. Measures gaussianity of data. ### stack all deltas on top of one another and measure dispersion that way + distribution.
		"""
		if hist:	
		
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
			plt.savefig(self.savefile + 'Run1a_deltas_distribution' + self.extension, dpi = self.dpi)
			
			if 'show' in list(kwargs.keys()) and kwargs['show']==True:
				plt.show()
			
			if 'params' in kwargs:
				params = kwargs['params']
				pickle.dump(ax, open(os.getcwd() + '/oo_analysis/figures/pickles/deltas_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
			else:
				pickle.dump(ax, open(os.getcwd() + "/oo_analysis/figures/pickles/deltas.pickle", 'wb'))
		else:
			deltas = self.grand_spectra_group['power_deviation'][...]
			
			ax = plt.subplot(111)
			plt.plot(deltas)
			plt.ylabel("power deviation (W)")
			plt.xlabel("Bins")
			
			plt.savefig(self.savefile + "Run1a_deltas" + self.extension, dpi = self.dpi)
			
			if 'show' in list(kwargs.keys()) and kwargs['show']==True:
				plt.show()
			
			if 'params' in kwargs:
				params = kwargs['params']
				pickle.dump(ax, open(os.getcwd() + '/oo_analysis/figures/pickles/deltas_distribution_(' + params['filter'] + "_" + params['pec_vel'] + "_" + params['signal'] + ").pickle", 'wb'))
			else:
				pickle.dump(ax, open(os.getcwd() + "/oo_analysis/figures/pickles/deltas_distribution.pickle", 'wb'))
		
		plt.clf()		
	def axion_model_compare(self, **kwargs):
		"""
			
		"""
	
		plt.clf()	
	def temp_v_freq(self, **kwargs):
		"""
		Plot the temperature over frequency
		
		
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
	def bg_sizes(self, **kwargs):
		"""
		Plot the background size, defined as the relative size of the background to the noise size, versus the mode frequency
		"""
		frequencies = []
		if 'mode_frequencies' not in kwargs.keys():
			for key in self.digitizer_group:
				attrs = self.digitizer_group[key].attrs
				if 'alog_timestamp' in attrs and 'cut' in attrs and not attrs['cut']:
					frequencies.append(attrs['mode_frequency'])		
		else:
			frequencies = kwargs['mode_frequencies']
			
		assert 'background_sizes' in list(kwargs.keys()) #Make sure background size is passed to kwargs
		codomain = kwargs['background_sizes']
		
		#Sort frequencies
		sort_mask = np.argsort(frequencies)
		master_mask = self.RF_interference_mask
		frequencies = np.array(frequencies)[sort_mask]
		bg_sizes = np.array(codomain)[sort_mask]
		
		fig = plt.figure()
		plt.xlabel("Mode Frequencies (MHz)")
		plt.ylabel("Relative Background Size")
		plt.plot(frequencies, bg_sizes)
		
		plt.savefig(self.savefile + "bg_sizes" + self.extension, dpi = self.dpi)
		
		plt.clf()			
	def synth_injection(self, frequency, **kwargs):
		"""
		Plot the background-removed scans that contain a synthetically injected axion signal. 
		Also plot the same portion of the axion fit significance
		"""
		
		afreq = str(frequency)
		self.delta_group = self.grand_spectra_group['deltas']
		injected_datasets_dict = {key: self.digitizer_group[key] for key in self.delta_group if 'Synthetic axion injected' in str(self.digitizer_group[key].attrs['notes']) and afreq in str(self.digitizer_group[key].attrs['notes'])}
		
		self.injected_datasets = list(injected_datasets_dict.values())[::2]
		injected_keys = list(injected_datasets_dict.keys())[::2]
		generator = range(len(injected_keys))
		
		fstarts = [float(i.attrs['start_frequency']) for i in self.injected_datasets]
		fstops = [float(i.attrs['stop_frequency']) for i in self.injected_datasets]
		ress = [float(i.attrs['frequency_resolution']) for i in self.injected_datasets]
		
		
		domains = {str(injected_keys[i]): np.arange(fstarts[i], fstops[i], ress[i]) for i in generator}
		offsets = {injected_keys[i]: 0.05*i for i in generator}
		
		plt.figure()
		plt.subplot(211)
		plt.yticks([])
		plots = [plt.plot(domains[i], self.delta_group[i][...] + offsets[i], color='black', alpha = 0.7, linewidth=0.5) for i in injected_keys]
		plt.xlabel("Frequency (MHz)")
		plt.subplot(212)
		self.data = self.grand_spectra_group['optimal_weight_sum'][...]/self.grand_spectra_group['model_excess_sqrd'][...]/self.grand_spectra_group['axion_fit_uncertainty'][...]
		mask = np.isfinite(self.data)
		mask1 = (self.data>0)
		self.master_mask = mask&self.RF_interference_mask&mask1
		reduced_data = self.data[self.master_mask]
		m = min([min(i) for i in list(domains.values())])
		M = max([max(i) for i in list(domains.values())])
		domain = self.axion_frequencies_MHz[self.master_mask]
		inx_m = np.argmin(np.abs(domain-m))
		inx_M = np.argmin(np.abs(domain-M))
		
		red_domain = domain[inx_m:inx_M]
		red_data = reduced_data[inx_m:inx_M]
		
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel("Significance")	
		plt.tight_layout()
		plt.plot(red_domain, red_data, label='Fit Significance')
			
		
		plt.savefig(self.savefile + "synth_injection({0})".format(afreq) + self.extension, dpi = self.dpi)
		
		if 'show' in list(kwargs.keys()) and kwargs['show']==True:
			plt.show()
		plt.clf()			
	def modulation_effect(self, **kwargs):
		"""
		Plot the variation of the two axion signals (N-Body and Standard Halo Model) as they change annually.
		"""
		h=4.135*10**(-15)
		c=2.998*10**5
		RME = 10**(-6) #RME of axion in eV
		sig = 160 #velocity dispersion of axion
		freq_o = RME/(h)
		f = np.linspace(freq_o, freq_o+500, 1000)
		delta_freq = f - freq_o
		energy = f*h     #Energy of equivalent photon
		kin_energy = energy-RME     #kinetic energy of axion
		vel_axion = c*np.sqrt(2*(kin_energy/RME))



		solar_vel=232.0 #km/s
		orbital_vel= 30.0 #km/s
		rotation_vel=0.46 #km/s

		vel_solar = solar_vel
		vel_min = solar_vel - orbital_vel - rotation_vel   #The lab velocity is antiparallel with the solar velocity.
		vel_max = solar_vel + orbital_vel + rotation_vel  #The lab velocity is parallel to the solar velocity

		x_coefficient_solar = (2*h*c**2)/(np.sqrt(2*np.pi*sig**2)*RME*vel_solar)
		x_coefficient_min = (2*h*c**2)/(np.sqrt(2*np.pi*sig**2)*RME*vel_min)
		x_coefficient_max = (2*h*c**2)/(np.sqrt(2*np.pi*sig**2)*RME*vel_max)

		x_coefficient=np.array([x_coefficient_solar, x_coefficient_min, x_coefficient_max])

		y_solar = np.sinh((vel_axion*vel_solar)/(sig**2))
		y_min = np.sinh(vel_axion*vel_min/(sig**2))
		y_max = np.sinh(vel_axion*vel_max/(sig**2))

		y_factor=np.array([y_solar, y_min, y_max])

		z_solar = np.exp(-(vel_axion**2 + vel_solar**2)/(2*sig**2))
		z_min = np.exp(-(vel_axion**2 + vel_min**2)/(2*sig**2))
		z_max = np.exp(-(vel_axion**2 + vel_max**2)/(2*sig**2))

		z_factor = np.array([z_solar, z_min, z_max])

		SHM_function_solar = x_coefficient[0] * y_factor[0] * z_factor[0]
		SHM_function_min = x_coefficient[1] * y_factor[1] * z_factor[1]
		SHM_function_max = x_coefficient[2] * y_factor[2] * z_factor[2]
		
		
		alphaexp = 0.36
		beta = 1.39			#cite values of Lentz et al. 2017
		T = 4.7*10**-7
		m = RME/(c**2) # Rest mass of axion in ev/c^2
		v_a = solar_vel #speed of sun in km/s
		h = 4.135*10**-15	#ev*s
		vo = freq_o #rest mass freq in MHz, about 241.83 MHz
		
		dv_min = -30.46 #km/s
		dv_max = +30.46 # km/s
		BoostMin = 1-(dv_min*v_a)/(v_a**2)
		BoostMax = 1-(dv_max*v_a)/(v_a**2)

		
		gamma = 1.012927949676458 #value of gamma function with input (1+alpha)/beta
		Cnum = beta
		Cden = ((1/(T*vo))**beta)**(-(1.0+alphaexp)/beta)*(1.0/(T*vo))**(alphaexp)*gamma
		Cden_min = ((BoostMin/(T*vo))**beta)**(-(1.0+alphaexp)/beta)*(BoostMin/(T*vo))**(alphaexp)*gamma
		Cden_max = ((BoostMax/(T*vo))**beta)**(-(1.0+alphaexp)/beta)*(BoostMax/(T*vo))**(alphaexp)*gamma
		C = Cnum/Cden
		C_min = Cnum/Cden_min
		C_max = Cnum/Cden_max
		
		MAXWELLFORM = C*((h*(f-vo)/(RME*T))**alphaexp)*np.exp(-(h*(f-vo)/(RME*T))**beta)
		MAXWELLFORM_min = C_min*((h*((f-vo)*BoostMin)/(RME*T))**alphaexp)*np.exp(-(h*((f-vo)*BoostMin)/(RME*T))**beta)
		MAXWELLFORM_max = C_max*((h*((f-vo)*BoostMax)/(RME*T))**alphaexp)*np.exp(-(h*((f-vo)*BoostMax)/(RME*T))**beta)
		

		plt.figure()
		plt.xlabel('Frequency (MHz)')
		plt.ylabel('Probability Density')
		plt.title("SHM and NBODY lineshape")
		
		plt.xlim((np.min(f)-5)*10**(-6), np.max(f)*10**(-6))
		plt.ylim(0, .008)
		plt.plot(f*10**(-6), MAXWELLFORM, linestyle='-.', linewidth = self.linewidth, alpha = 1, color='green', label='N-Body Solar Average')
		plt.plot(f*10**(-6),SHM_function_solar, linestyle=':', linewidth = self.linewidth, alpha = 1, color='black', label='SHM Solar Average')
		plt.fill_between(f*10**(-6), MAXWELLFORM, MAXWELLFORM_max, facecolor = 'blue', alpha = 0.35)
		plt.fill_between(f*10**(-6), MAXWELLFORM, MAXWELLFORM_min, facecolor = 'red', alpha = 0.35)
		plt.fill_between(f*10**(-6), SHM_function_solar, SHM_function_min, facecolor = 'red', alpha = 0.35)
		plt.fill_between(f*10**(-6), SHM_function_solar, SHM_function_max, facecolor = 'blue', alpha = 0.35)
		plt.ticklabel_format(useOffset=False)
		plt.yticks(np.arange(0, .008, .001))
		plt.xticks(np.arange(vo*10**(-6), (vo+500)*10**(-6), 200*10**(-6)), rotation = 20)
		plt.ticklabel_format(style='plain')
		
		plt.legend(loc=1, frameon=True, framealpha=1.0)
		plt.tight_layout()
		
		plt.savefig(self.savefile + "modulation_effect" + self.extension, dpi = self.dpi)
		plt.show()	
	def sensitivity_combined(self, **kwargs):
		"""
		Combined coupling sensitivity, Power Sensitivity, and dark matter sensitivity using the SHM and the NBody signals
		"""
		import pickle
		import pandas as pd
		nbody_files = ['oo_analysis/figures/pickles/Coupling_sensitivity_(RCHPF_earth rotation_axionDM_w_baryons).pickle', 'oo_analysis/figures/pickles/Power_sensitivity_(RCHPF_earth rotation_axionDM_w_baryons).pickle', 'oo_analysis/figures/pickles/DM_sensitivity_(RCHPF_earth rotation_axionDM_w_baryons).pickle']
		
		shm_files = ['oo_analysis/figures/pickles/Coupling_sensitivity_(RCHPF_earth rotation_SHM).pickle', 'oo_analysis/figures/pickles/Power_sensitivity_(RCHPF_earth rotation_SHM).pickle', 'oo_analysis/figures/pickles/DM_sensitivity_(RCHPF_earth rotation_SHM).pickle']
			
		fig = plt.figure()
		nbody_pickles = [pickle.load(open(i, 'rb')) for i in nbody_files]
		
		shm_pickles = [pickle.load(open(i, 'rb')) for i in shm_files]
		
		get_ydata = lambda obj: np.array(obj.get_lines()[0].get_ydata())
		get_xdata = lambda obj: np.array(obj.get_lines()[0].get_xdata())
		
		nbody_ydata = [get_ydata(i) for i in nbody_pickles]
		nbody_xdata = [get_xdata(i) for i in nbody_pickles]
		
		shm_ydata = [get_ydata(i) for i in shm_pickles]
		shm_xdata = [get_xdata(i) for i in shm_pickles]
		
		
		
		#Coupling_sensitivity
		plt.clf()
		plt.cla()
		plt.show()
		plt.close()
		
		
		plt.figure(figsize = (8.1,6))
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel(r"Coupling Sensitivity ($g_{a \gamma \gamma}/g_{a \gamma \gamma, DFSZ}$)")
		
		plt.locator_params(axis='x', nbins=10)
		plt.locator_params(axis='y', nbins=1)
		plt.tick_params(which='minor', length=2.5)
		plt.tick_params(which='major', length=5.0)
		plt.xticks(rotation=25)
		plt.ticklabel_format(useOffset=False)
		plt.yscale('log')
		plt.plot(nbody_xdata[0], nbody_ydata[0], label='N-Body', color='blue', ls='steps')
		plt.plot(shm_xdata[0], shm_ydata[0], label='SHM', color='green', ls='steps')
		plt.axhline(y=1, linestyle='--', color='black',  label='DFSZ Coupling')
		plt.axhline(y=0.97/0.35, linestyle = ':', color='black', label='KSVS Coupling')
		plt.legend(loc=1, frameon=True, framealpha=1.0)
		plt.savefig(self.savefile + "combined_couplings" + self.extension, dpi = self.dpi)
		if self.save_to_paper:
			plt.savefig(self.paper_dir +  "combined_couplings" + self.extension, dpi = self.dpi) #Save image
		if 'show' in kwargs.keys() and kwargs['show']:
			plt.show()
		
		plt.clf()
		plt.cla()
		
		
		
		#Power Sensitivity
		plt.figure(figsize = (8.1,6))

		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel(r"Power Sensitivity ($P/P_{DFSZ}$)")
			
		plt.locator_params(axis='x', nbins=10)
		plt.locator_params(axis='y', nbins=1)
		plt.tick_params(which='minor', length=2.5)
		plt.tick_params(which='major', length=5.0)
		plt.xticks(rotation=25)
		plt.ticklabel_format(useOffset=False)
		plt.yscale('log')
		plt.plot(nbody_xdata[1], nbody_ydata[1], label='N-Body', color='blue', ls='steps')
		plt.plot(shm_xdata[1], shm_ydata[1], label='SHM', color='green', ls='steps')
		plt.axhline(y=1, linestyle='--', color='black',  label='DFSZ Power')
		plt.axhline(y=(0.97/0.36)**2, linestyle = ':', color='black', label='KSVZ Power')
		plt.legend(loc=1, frameon=True, framealpha=1.0)
		plt.savefig(self.savefile + "combined_power" + self.extension, dpi = self.dpi)
		if self.save_to_paper:
			plt.savefig(self.paper_dir +  "combined_power" + self.extension, dpi = self.dpi) #Save image
		if 'show' in kwargs.keys() and kwargs['show']:
			plt.show()
		plt.clf()
		
		
		#Dark Matter Sensitivity
		plt.figure(figsize = (8.1,6))
			
		plt.xlabel("Axion Frequencies (MHz)")
		plt.ylabel(r"Dark Matter Sensitivity ($DFSZ/GeV/cc$)")
		
		
		plt.ylim(1.2*10**(-1), 1)
		plt.locator_params(axis='x', nbins=10)
		plt.locator_params(axis='y', nbins=1)
		plt.tick_params(which='minor', length=2.5)
		plt.tick_params(which='major', length=5.0)
		plt.xticks(rotation=25)
		plt.ticklabel_format(useOffset=False)
		plt.yscale('log')
		plt.plot(nbody_xdata[2], nbody_ydata[2], label='N-Body', color='blue', ls='steps')
		plt.plot(shm_xdata[2], shm_ydata[2], label='SHM', color='green', ls='steps')
		#plt.axhline(y=0, linestyle='--', color='black',  label='DFSZ Power')
		plt.axhline(y=0.45, linestyle = ':', color='black', label='Expected Halo Density')
		plt.legend(loc=1, frameon=True, framealpha=1.0)
		plt.savefig(self.savefile + "combined_DM" + self.extension, dpi = self.dpi)
		if self.save_to_paper:
			plt.savefig(self.paper_dir +  "combined_DM" + self.extension, dpi = self.dpi) #Save image
		if 'show' in kwargs.keys() and kwargs['show']:
			plt.show()
		plt.clf()
	def parsed_artificial_scan(self, **kwargs):
		"""
		Generate an artificial scan and parse out the relevant components (Background, noise, axion signal). Then plot the fourier transform in [power]/[s].
		Next plot the same artificial scan but extended using the method used in the Reciprocated clone HPF. THen plot its fourier transform in [power]/[s]
		"""
		signal_color = 'yellow'
		axion_color = 'blue'
		background_color = 'brown'
		noise_color='green'

		from oo_analysis.filters.RCHPF import generate_signal, gen_recip_copy
		n = 5
		
		sig = generate_signal() #Dictionary containing all the components of the artificial scan (Total, Background, noise, axion signal)
		noise = sig['noise']
		signal = sig['Combined signal']
		axion = sig['axion signal']
		back = sig['backgrounds']
		
		autocorr = lambda arr: np.correlate(arr-np.mean(arr), arr-np.mean(arr), 'full')[len(arr)-1:]
		fourier = lambda arr: np.abs(np.fft.fftn(arr**2))[:int(len(arr)/2)]
		noise_corr = fourier(noise)
		axion_corr = fourier(axion)
		back_corr = fourier(back)/fourier(back)[0]

		#Plot signal and its fourier transform
		fig, ax = plt.subplots(2,2)
		ax[0,0].plot(signal, color=signal_color, label='Signal'); 
		ax[0,0].plot(back, color=background_color, label='Background'); 
		ax[0,0].plot(noise, color=noise_color, label='noise'); 
		ax[0,0].plot(axion, color=axion_color, label='Axion signal'); 
		ax[0,0].legend(loc='upper right'); 
		ax[0,0].set_xlabel("Bins"); 
		ax[0,0].set_ylabel("Amplitude"); 


		ax[0,1].plot(back_corr*10, color=background_color, label='Background'); 
		ax[0,1].plot(noise_corr*50, color=noise_color, label='noise'); 
		ax[0,1].plot(axion_corr*10**(46), color=axion_color, label='Axion signal'); 
		ax[0,1].legend(loc='upper right'); 
		ax[0,1].set_xlabel("Bins"); 
		ax[0,1].set_ylabel("Amplitude"); 
		#plt.ylim(-0.01, 1); 
		#plt.xlim(-0.25, 256/2); 
		ax[0,1].set_yscale("log") 

		#Now extend the artificial scan by cloning, flipping, and reciprocating. PLots it and its fourier transform
		ext_signal = gen_recip_copy(signal, n); 
		ext_noise = gen_recip_copy(noise+1, n)-1; 
		ext_axion = gen_recip_copy(axion+1, n)-1; 
		ext_back = gen_recip_copy(back, n); 
		ext_fourier_axion = fourier(ext_axion); 
		ext_fourier_back = fourier(ext_back)/fourier(ext_back)[0]; 
		ext_fourier_noise = fourier(ext_noise); 


		ax[1,0].plot(ext_signal, color=signal_color, label='Signal'); 
		ax[1,0].plot(ext_noise, color=noise_color, label='Noise'); 
		ax[1,0].plot(ext_axion, color=axion_color, label='Axion Signal'); 
		ax[1,0].plot(ext_back, color=background_color, label='Background'); 

		ax[1,0].legend(loc='upper right');  
		ax[1,0].set_xlabel("Bins"); 
		ax[1,0].set_ylabel("Amplitude");



		ax[1,1].plot(ext_fourier_axion*10**(22), color=axion_color, label='Axion Signal'); 
		ax[1,1].plot(ext_fourier_noise*50, color=noise_color, label='Noise', alpha = 0.5); 
		ax[1,1].plot(ext_fourier_back*10, color=background_color, label='Background', alpha = 0.5); 
		ax[1,1].legend(loc='upper right'); 
		ax[1,1].set_xlabel('Bins'); 
		ax[1,1].set_ylabel("Amplitude");
		#plt.xlim(0, len(ext_fourier_noise)/2); 
		#plt.ylim(-0.01,1); 
		ax[1,1].set_yscale("log")


		plt.tight_layout()
		plt.savefig(os.getcwd() + "/oo_analysis/figures/parsed_extended_domain_ffts" + self.extension, dpi = self.dpi)
		
		if 'show' in kwargs.keys() and kwargs['show']:
			plt.show()
		plt.clf()
	def extended_fourier_domain(self, **kwargs):
		"""
		Makes a plot of the extended fourier domain
		"""
		from oo_analysis.filters.RCHPF import generate_signal, gen_recip_copy
		
		n = 1
		sig = generate_signal()
		back = sig['backgrounds']		
		ext_back = gen_recip_copy(back, n); 
		
		plt.title("Extended Fourier Domain")
		plt.xlabel("bins")
		plt.ylabel("Amplitude")
		plt.plot(ext_back)
		[plt.axvline((i+1)*256, color='black', linewidth=1, linestyle=':') for i in range(3)]
		
		upper_pos = (max(back)-back[0])/2 + back[0]
		lower_pos = -(back[0]-min(ext_back))/2 + back[0]
		
		plt.text(33, upper_pos, "Reciprocated \n      clone")
		plt.text(280, lower_pos, "Original signal")
		plt.text(545, upper_pos, 'Reciprocated \n      clone')
		plt.text(801, lower_pos, "Cloned signal")
		plt.xlim(0, len(ext_back)+1)
		
		plt.savefig(os.getcwd() + "/extended_fourier_domain" + self.extension, dpi = self.dpi)
		if 'show' in kwargs.keys() and kwargs['show']:
			plt.show()
		plt.clf()
	def change_in_coup_sens(self, **kwargs):
		"""
		Plots the change in coupling sensitivity between using a dynamical lineshape and a static lineshape
		"""
		
		import pickle
		import pandas as pd
		
		earth_rot_files = os.getcwd() + '/oo_analysis/figures/pickles/Coupling_sensitivity_(RCHPF_earth rotation_SHM).pickle'
		sol_files = os.getcwd() + '/oo_analysis/figures/pickles/Coupling_sensitivity_(RCHPF_solar_SHM).pickle'
		
		rot_pickles = pickle.load(open(earth_rot_files, 'rb'))
		
		sol_pickles = pickle.load(open(sol_files, 'rb'))
		
		get_ydata = lambda obj: np.array(obj.get_lines()[0].get_ydata())
		get_xdata = lambda obj: np.array(obj.get_lines()[0].get_xdata())
		
		rot_ydata = get_ydata(rot_pickles)
		rot_xdata = get_xdata(rot_pickles)
		
		sol_ydata = get_ydata(sol_pickles)
		sol_xdata = get_xdata(sol_pickles)
		
		difference = rot_ydata - sol_ydata
		domain = rot_xdata
		
		fig, ax = plt.subplots(1)
		ax[0,0].set_xlabel("Axion Frequencies (MHz)")
		ax[0,0].set_ylabel(r"Coupling Sensitivity Difference ($g_{a \gamma \gamma}/g_{a \gamma \gamma, DFSZ}$)")
		ax[0,0].plot(domain, difference)
		
		plt.savefig(os.getcwd() + "/oo_analysis/figures/coupling_sens_difference" + self.extension, dpi = self.dpi)
		
		if self.save_to_paper:
			plt.savefig(self.paper_dir +  "Rotation_minus_solar_coupling_sensitivity" + self.extension, dpi = self.dpi) #Save image
		if 'show' in kwargs.keys() and kwargs['show']:
			plt.show()
		plt.clf()





		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		