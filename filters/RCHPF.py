##System management
import os
import sys
sys.path.append("..")
##


##Packages for arithmetic, constants, and general pythonic containers
import numpy
import astropy.constants as C
import astropy.units as u
import math
from math import gamma
import copy
##

##Background subtraction specific packages
import filters.backsub_filters_lib
from filters.backsub_filters_lib import SG, poly_fit #Savitzky Golay Filter as 6 order polynomial fit
from scipy.stats import variation #Coefficient of Variation
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
#import argparse
##

## Import packages for parallel-processing
import pyopencl as cl 
import reikna.cluda as cluda
from reikna.fft import FFT as fftcl
##

## Import packages for fourier transforms
from toolbox.DFT import DFT, IDFT
##

##Meta Analysis
import time
##

##Plotting software
from toolbox.plot_dataset import plotter
import matplotlib.pyplot as plt
import pandas as pd

##


os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1" #Output pyopencl compiler messages
	
def reciprocated_clone_hpf(data, npairs, testing=False, testing_subdir = "", scan_number="",iter = '', **meta):
	"""
	Parameters:
				data: can be array of data or a digitizer number
				npairs: number of clones to make when generating large fourier domain
				(removed)ch: --OBSOLETE-- specify channel data was taken on.
				testing: Used for testing... Plots major calculations in the script of meta analysis if True
				scan_number: used for testing. Specifies scan number in title when saving plots
	"""
	ORIG = copy.deepcopy(data)
	avail_types = [numpy.ndarray, numpy.arange(1), 1, 1.]
	for i in avail_types:
		if isinstance(data, i):
			break
		else:
			return "Error: invalid data type"
	errors = {'maxed_filter_size':False}
	
	#if (numpy.mean(data)-data[0])>numpy.std(data): 
	ORIG[0]=numpy.mean([ORIG[1], ORIG[2], ORIG[3]]) # Some spectra have odd behavior at beginning of scan, i.e. a single downward spike at the beginning position. I just set a default value
	
	while len(ORIG)<256: #For some reason, some datasets dont have 2**8 elements.
		ORIG = numpy.append(ORIG[-1], ORIG)
	
	reflecting_time_start = time.time()
	sz = len(ORIG)
	n = npairs # Specifies number of concatenated datasets
	reflecting_time_stop = time.time()
	
	#Generate statistics about the data set and define variables and functions
	
		#reflects the dataset about the center and operates on it to sequence datasets without discontinuities 

	calculating_tophat_size_start = time.time()
	reciprocating_array_start = time.time()

	RECIP, lin_fit_beg, lin_fit_end = gen_recip_copy(ORIG)
	reciprocating_array_stop = time.time()
	
	#Construct large array consisting of reflected and nonreflected data_array
	generating_large_array_start = time.time()
	appended_ORIG = [*RECIP, *ORIG]*(n+1)
	szA = len(appended_ORIG) # Find length of large array
	generating_large_array_stop = time.time()

	sigma = calc_filter_size(ORIG, appended_ORIG)
	if sigma>2**8:
		print("Ceiling of reciprocated clone hpf reached (sigma={0:0.5f})".format(sigma))
		sigma=2**8
		errors['maxed_filter_size'] = True
		
	calculating_tophat_size_stop = time.time()
		

	"""
	tophat = [0]*szA # Generate tophat function to pick out lower frequency structure
	for i in range(int(sigma)):
		tophat[i]=1
		tophat[szA-i-1] = 1
	"""
	
	generating_tophat_start = time.time()
	#experimental tophat
	tophat_func = lambda x, sigma, size: 1/(1 + (x/sigma)**12) + 1/(1 + ((x-size)/sigma)**12)
	tophat = tophat_func(numpy.arange(szA), sigma, szA)
	generating_tophat_stop = time.time()
	
	
	# Fourier transform large array
	fft_and_highpassfilter_start = time.time()
	fft_appended_ORIG = DFT(appended_ORIG) 
	reduced_fft_appended_ORIG = tophat*fft_appended_ORIG # Isolate lower frequencies
	fft_and_highpassfilter_stop = time.time()
	
	
	# Inverse Fourier Transform
	ifft_start = time.time()
	BS_fft_appended_ORIG = IDFT(reduced_fft_appended_ORIG)
	ifft_stop = time.time()
	
	#pick out the original signal
	
	picking_og_signal_start = time.time()
	pickORIG = [0]*len(ORIG)
	if n%2==1: #calculate which part of the array to pick out. if n is odd, pick out the original scan to the left of the center of the array. if n is even, pick out the middle
		l=n
	elif n%2==0:
		l=n+1
		
	for i in range(len(ORIG)):
		pickORIG[i] = BS_fft_appended_ORIG[(l*len(ORIG)+i)]
	picking_og_signal_stop = time.time()
	
	dividing_structure_start = time.time()
	ORIG = numpy.array(ORIG) #Convert back into an array
	pickORIG =  numpy.array(pickORIG)
	filtereddata = ORIG/pickORIG # Divide out low freq. structure from data.
	dividing_structure_stop = time.time()
	
	if testing:
		import datetime as dt
		date = "(" + dt.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + ")"
		starting_string=""
		if scan_number!="":
			starting_string = '('+str(scan_number)+')_'
		sdir = "D:/Users/shaun/Documents/Coding Software/Python/Scripts/New-Analysis-Scheme/oo_analysis/figures/Background_sub_testing/rchpf_testing/"+ str(testing_subdir)+starting_string
		format = '.pdf'
		n="("+str(iter)+")"
		ylim = (min(ORIG), max(ORIG)*(1+1/24))
		plotter(ORIG.real,title = n+'data', ylimits=ylim, savedir = sdir+n+'data' + date +format); 
		plotter(ORIG[::-1].real, title = n+'refl_data', ylimits=ylim, savedir = sdir+n+'refl_data' + date+format); 
		plotter(RECIP.real, title=n+'recip_data', savedir = sdir+n+'recip_datas' + date+format); 
		plotter(appended_ORIG, title=n+'large_append_array', savedir = sdir+n+'large_append_array' + date+format); 
		plotter(fft_appended_ORIG.real, title = n+'fft-ed large append array' + "With sigma = {0}".format(sigma), savedir = sdir+n+'fft-ed_large_appended_array' + date+format); 
		plotter(BS_fft_appended_ORIG.real, title = n+'inverse_fft-ed_append_array', savedir = sdir+n+'inverse_fft-ed_append array' + date+format);
		plotter(pickORIG.real, title = n+'recovered data from ifft-ed array', savedir = sdir+n+'recovered_data_from_ifft-ed_array' + date+format)
		plotter(filtereddata.real, title = n+'filtered data ' + 'sigma={0:0.3f}'.format(sigma), ylimits = (min(filtereddata), max(filtereddata)*(1+1/24)), savedir = sdir+n+'filtered_data' + date+format)
		plotter(tophat, title = n+'filter function' ,savedir = sdir + n+"filter_function" + date+format)
	
	#meta analysis
	if 'timeit' in meta.keys():
		if meta['timeit']:
			reflecting_time = reflecting_time_stop - reflecting_time_start
			calculating_tophat_size = calculating_tophat_size_stop - calculating_tophat_size_start
			reciprocating_array = reciprocating_array_stop - reciprocating_array_start
			generating_large_array = generating_large_array_stop - generating_large_array_start
			generating_tophat = generating_tophat_stop - generating_tophat_start
			fft_and_highpassfilter = fft_and_highpassfilter_stop - fft_and_highpassfilter_start
			ifft =ifft_stop- ifft_start
			picking_og_signal = picking_og_signal_stop - picking_og_signal_start
			dividing_structure = dividing_structure_stop - dividing_structure_start
			meta['reflecting_time'].append(reflecting_time)
			meta['calculating_tophat_size'].append(calculating_tophat_size)
			meta['reciprocating_array'].append(reciprocating_array)
			meta['generating_large_array'].append(generating_large_array)
			meta['generating_tophat'].append(generating_tophat)
			meta['fft_and_highpassfilter'].append(fft_and_highpassfilter)
			meta['ifft'].append(ifft)
			meta['picking_og_signal'].append(picking_og_signal)
			meta['dividing_structure'].append(dividing_structure)
	
	return {"filtereddata":filtereddata.real, "sigma": sigma, "number of clones": npairs, "meta": meta}, errors

def gen_recip_copy(arr):
	num_fit_points = 15
	linear_func = lambda a,b,arr: a*arr + b
	domain = numpy.arange(num_fit_points)
	
	popt_beginning = numpy.polyfit(domain, arr[0:num_fit_points], 1)
	lin_fit_beginning_clone = linear_func(popt_beginning[0], popt_beginning[1], domain)
	
	popt_end = numpy.polyfit(domain, arr[len(arr)-num_fit_points:], 1)
	lin_fit_end_clone = linear_func(popt_end[0], popt_end[1], domain)
	
	#RECIP = ((arr[-1]**2 - arr[0]**2)*(numpy.cos(numpy.arange(len(arr))*numpy.pi/(2*numpy.size(arr[::-1])))**2) + arr[0]**2)/arr[::-1]
	RECIP = ((lin_fit_end_clone[-1]**2 - lin_fit_beginning_clone[0]**2)*(numpy.cos(numpy.arange(len(arr))*numpy.pi/(2*len(arr)))**2) + lin_fit_beginning_clone[0]**2)/arr[::-1]

	
	
	return RECIP, lin_fit_beginning_clone, lin_fit_end_clone

def calc_filter_size(arr, recip_arr):
	
	### Currently testing different calculations of the filter size
	delta = lambda arr: (numpy.max(arr)-numpy.min(arr))
	arg_delta = lambda arr: numpy.abs((numpy.argmax(arr)-numpy.argmin(arr)))

	sigma_func = lambda arr: delta(arr)/numpy.std(arr[:10])
	
	ORIG_corr = numpy.correlate(recip_arr-numpy.mean(recip_arr), recip_arr-numpy.mean(recip_arr), 'full')
	
	
	#plt.plot(ORIG_corr); plt.show(); plt.plot(recip_arr); plt.show();

	#Find width of structure using the autocorrelation
	peaks,_ = find_peaks(ORIG_corr, height=numpy.max(ORIG_corr)/2, distance=50)
	primary_peak = peaks[numpy.where(ORIG_corr[peaks]==numpy.max(ORIG_corr))[0]][0]
	peaks_rem = numpy.delete(peaks, numpy.where(peaks==primary_peak))
	secondary_peak = peaks_rem[numpy.max(numpy.where(ORIG_corr[peaks_rem]==numpy.max(ORIG_corr[peaks_rem]))[0])]
	width = secondary_peak-primary_peak
	width_factor = 1/numpy.tanh(width/(512/2))
	
	sigma = sigma_func(arr)*width_factor
	return sigma














####################################
#Scripts for analyzing different background subtraction methods
####################################


#############
#savedir = "C:/Users/drums/Documents/ADMX/Papers/Own Papers/Plots and Graphics/Bsub/MyBSub/"
#############


# data and axion params 
nsamp = 1
c = C.c.to(u.km/u.s).value # in km/s
h = C.h.to(u.eV*u.s).value # in eV*s
k = C.k_B.value #J/K
ma = 1.e-6 #eV
rest_energy = ma
nbins = 2**8 # Length of scan
df = 95.4 #Bin Size

#Construct axion frequency array
start = rest_energy/h - nbins*df/2
stop = start + nbins*df
bins = [start + df*i for i in range(nbins)]

def ADMX_signal_model(freq,exp_pow,poly_pow,temp):
    if freq<=rest_energy/h:
        return 0.0
    else:
        pass
    sig = temp*rest_energy/h
    nuosig = (freq-rest_energy/h)/sig
    gamma = math.gamma((1.0+poly_pow)/exp_pow)
    Cnum = exp_pow
    Cden = ((1.0/(sig))**exp_pow)**(-(1.0+poly_pow)/exp_pow)*(1.0/sig)**(poly_pow)*gamma
    C = Cnum/Cden
    dist = C*(nuosig)**poly_pow*numpy.exp(-(nuosig)**exp_pow)
    return dist
signal_strength = 0.15/3 # Level to be recoverable with ~20 scans at SNR=3
params = [1.4,0.37,4.7e-7]

#####Parameters specific to background subtraction testing
sig_atten_array = [["sg"],["bs"],["poly"],["noise"]]
dispersion_array = [["sg"],["bs"],["poly"],["noise"]]
power_spec_array = [["sg"],["bs"],["poly"],["noise"]]
axion_signal_array = []
autocorr_array = [[],[],[],[]]
savedir = "D:/Users/shaun/Documents/Coding Software/Python/Scripts/New-Analysis-Scheme/oo_analysis/figures/Background_sub_testing/"

def generate_signal(with_signal=True, noise_only=False, **kwargs):

	level = 10**(-9)
	signal_strength = 0.05/3 # Level to be recoverable with ~20 scans at SNR=3
	params = [1.4,0.37,4.7e-7]
	order = 3
	#random polynomials
	coeffs = (2*numpy.random.random((order,))-1.0)
	rand_poly = lambda x: sum([coeff*x**i for i,coeff in enumerate(coeffs)])
	#random hot rod
	Q = numpy.random.uniform(low=10000, high=60000)
	coeff = numpy.random.uniform(low=0.0, high=1)
	rand_hot_rod = lambda x: coeff/(1.0+(x-rest_energy/h)**2/(rest_energy/(h*Q))**2)
	
	polys=[]
	for i in range(nsamp):
		polys.append([0.1*rand_poly(1.0*x/nbins) for x in range(nbins)])
		
	ress = []
	for i in range(nsamp):
		ress.append([rand_hot_rod(x) for x in bins])
	background = numpy.ones((nsamp, nbins)) + numpy.array(polys) + numpy.array(ress)
	
	single_signal = numpy.array([signal_strength*df*ADMX_signal_model(x, params[0], params[1], params[2]) for x in bins])
	signals = numpy.array([single_signal for i in range(nsamp)])
	#Create noise test signals
	noise = 0.01*numpy.random.rand(nsamp, nbins)
	#combine everything
	if with_signal and not noise_only:
		Combined_Signal = numpy.array((noise+signals+1.0)*background)[0]*level
	elif noise_only and with_signal:
		Combined_Signal = numpy.array((noise+signals+1.0))[0]*level
	elif not with_signal and noise_only:
		Combined_Signal = numpy.array((noise+1.0))[0]*level
	noises = (noise+signals)
	toreturn = {"Combined signal":Combined_Signal, "noise":noises[0], "backgrounds":background[0], "axion signal":single_signal}
	return toreturn
    
def analysis(scan = [], noise_only=False,**kwargs):

	#Generate artificial data and push it through background subtraction functions
	generated_signal = generate_signal(noise_only = noise_only, **kwargs)
	sample_scan = generated_signal["Combined signal"]
	sample_noise = generated_signal["noise"]
	sample_backgrounds = generated_signal["backgrounds"]
	if len(scan)!=0:
		BS_res,_ = reciprocated_clone_hpf(scan, 3,  testing_subdir="BS_Analysis_runs/", **kwargs)
		BS_scan = BS_res["filtereddata"].real
		Poly_scan = poly_fit(scan) 
		SG_scan = SG(scan)
	else:
		BS_res,_ = reciprocated_clone_hpf(sample_scan, 3,  testing_subdir="BS_Analysis_runs/", **kwargs)
		BS_scan = BS_res["filtereddata"].real
		Poly_scan = poly_fit(sample_scan) 
		SG_scan = SG(sample_scan)
	nbins=len(sample_scan)
	signals = numpy.asarray([sample_scan, BS_scan, Poly_scan, SG_scan])
	

	
	

	power_spec_array[0].append((SG_scan-numpy.mean(SG_scan)))
	power_spec_array[1].append(BS_scan-numpy.mean(BS_scan))
	power_spec_array[2].append(Poly_scan-numpy.mean(Poly_scan))
	power_spec_array[3].append(sample_noise)
	axion_signal_array.append(generated_signal["axion signal"])
	
	#define plotting parameters and plot returned signals
	Mins=numpy.array([numpy.min([0.,numpy.min(sample_scan)]), numpy.min(sample_noise), numpy.min(BS_scan.real-1), numpy.min(Poly_scan-1.), numpy.min(SG_scan-1.)])*(1+(1/24))        #Min += -Min/24
	Maxs=numpy.array([numpy.max(sample_scan), numpy.max(sample_noise), numpy.max(BS_scan.real-1.), numpy.max(Poly_scan-1.), numpy.max(SG_scan-1.)])*(1+(1/12))       #Max += Max/12
	plotparams=["Bins", "Amplitude"]
	signal_to_plot = [sample_scan, sample_noise, BS_scan-1., Poly_scan-1., SG_scan-1.]
	titles = ["sample scan", "Sample scan noise characteristic", "RCHPF Background subtracted scan (deltas)", "six order polynomial fit (deltas)", "savitzky golay fit (deltas)"]
	filenames = ["sample_scan.pdf","sample_noise.pdf", "RCHPF_sample_scan.pdf","Poly_fit_scan.pdf","SG_fit_scan.pdf"]
	[plotter(signal_to_plot[i], title=titles[i], xlabel=plotparams[0], ylabel=plotparams[1],  ylimits = (Mins[i], Maxs[i]), savedir=savedir + filenames[i]) for i in numpy.arange(len(signal_to_plot))]


	#calculate and plot normalized autocorrelations
	pre_autocorrs = numpy.correlate(sample_scan, sample_scan, 'full')
	noise_autocorrs = numpy.correlate(sample_noise, sample_noise, 'full') 
	sg_autocorrs = numpy.correlate(SG_scan-1., SG_scan-1., 'full') #Notice correlations are computed on a mean-removed signal
	poly_autocorrs = numpy.correlate(Poly_scan-1., Poly_scan-1., 'full')
	BS_autocorrs = numpy.correlate(BS_scan-1., BS_scan-1., 'full')
	autocorrs = [pre_autocorrs, noise_autocorrs, sg_autocorrs, poly_autocorrs, BS_autocorrs]

	autocorr_array[0].append(sg_autocorrs[nbins:])
	autocorr_array[1].append(BS_autocorrs[nbins:])
	autocorr_array[2].append(poly_autocorrs[nbins:])
	autocorr_array[3].append(noise_autocorrs[nbins:])
	
	# calculate dispersion and mean of autocorrelations
	pre_autocorrs_stdev = numpy.std(pre_autocorrs[nbins:])
	noise_autocorrs_stdev = numpy.std(noise_autocorrs[nbins:])
	poly_autocorrs_stdev = numpy.std(poly_autocorrs[nbins:])
	BS_autocorrs_stdev = numpy.std(BS_autocorrs[nbins:])
	sg_autocorrs_stdev = numpy.std(sg_autocorrs[nbins:])

	pre_autocorrs_mean = numpy.mean(pre_autocorrs[nbins:])
	noise_autocorrs_mean = numpy.mean(noise_autocorrs[nbins:])
	sg_autocorrs_mean = numpy.mean(noise_autocorrs[nbins:])
	poly_autocorrs_mean = numpy.mean(noise_autocorrs[nbins:])
	BS_autocorrs_mean = numpy.mean(BS_autocorrs[nbins:])

	Min = numpy.min([numpy.min(sg_autocorrs[nbins:]), numpy.min(poly_autocorrs[nbins:]), numpy.min(BS_autocorrs[nbins:])])*(1-1/24)   #Min -= Min/24
	Max = numpy.max([numpy.max(sg_autocorrs[nbins:]), numpy.max(poly_autocorrs[nbins:]), numpy.max(BS_autocorrs[nbins:])])*(1+1/24)   #Max += Max/24
	length=len(pre_autocorrs[nbins:])

	autocorrs_to_plot = [pre_autocorrs[nbins:], noise_autocorrs[nbins:], sg_autocorrs[nbins:], poly_autocorrs[nbins:], BS_autocorrs[nbins:]]

	stdevs = [pre_autocorrs_stdev, noise_autocorrs_stdev, sg_autocorrs_stdev, poly_autocorrs_stdev, BS_autocorrs_stdev]
	
	dispersion_array[0].append(sg_autocorrs_stdev)
	dispersion_array[1].append(BS_autocorrs_stdev)
	dispersion_array[2].append(poly_autocorrs_stdev)
	dispersion_array[3].append(noise_autocorrs_stdev)
	
	means = [pre_autocorrs_mean, noise_autocorrs_mean, sg_autocorrs_mean, poly_autocorrs_mean, BS_autocorrs_mean]

	totals = [numpy.sum(numpy.abs(i)) for i in autocorrs_to_plot]

	autocorrs_titles = ["Signal Autocorrelation", "Signal Noise Autocorrelation", "Savitzky Golay Autocorrelation", "6-order Polynomial Fit Autocorrelation", "RCHPF Autocorrelation"]
	autocorrs_filenames = ["pre_autocorrs.pdf", "noise_autocorrs.pdf", "SG_autocorrs.pdf", "poly_autocorrs.pdf", "RCHPF_autocorrs.pdf"]
	autocorr_plotparams = ["1/bins", "Amplitude"]

	#[plotter(autocorrs_to_plot[i], title=autocorrs_titles[i], xlabel=autocorr_plotparams[0], ylabel=autocorr_plotparams[1], Min=np.min(autocorrs_to_plot[i]*(1-1/12)), Max=np.max(autocorrs_to_plot[i]*(1+1/12)), filename=autocorrs_filenames[i]) for i in np.arange(len(autocorrs_to_plot))]

	for i in numpy.arange(len(autocorrs_to_plot)):
		domain=numpy.arange(length)
		plt.ylabel(autocorr_plotparams[1])
		plt.xlabel(autocorr_plotparams[0])
		plt.title(autocorrs_titles[i])
		plt.plot(domain, autocorrs_to_plot[i])
		plt.text(x=150, y=max(autocorrs_to_plot[i])*(1-1/12), s="Mean = {0:0.5f} \n Stand. Dev = {1:0.5f} \n Totals = {2:0.5f}".format(means[i], stdevs[i], totals[i]))
		plt.tight_layout()
		plt.savefig(savedir+autocorrs_filenames[i], dpi=600)
		plt.clf()
	
	
	#signal-attenuation calculations and plots
	axion_signal = numpy.sum(generated_signal["axion signal"][128:133])
	SG_deltas = numpy.sum((SG_scan-numpy.mean(SG_scan))[128:133])
	BS_deltas = numpy.sum((BS_scan-numpy.mean(BS_scan))[128:133])
	POLY_deltas = numpy.sum((Poly_scan-numpy.mean(Poly_scan))[128:133])
	NOISE_deltas = numpy.sum((sample_noise-numpy.mean(sample_noise))[128:133])
	
	SG_sig_atten = SG_deltas/axion_signal
	BS_sig_atten = BS_deltas/axion_signal
	POLY_sig_atten = POLY_deltas/axion_signal
	NOISE_sig = NOISE_deltas/axion_signal
	sig_atten_array[0].append(SG_sig_atten)
	sig_atten_array[1].append(BS_sig_atten)
	sig_atten_array[2].append(POLY_sig_atten)
	sig_atten_array[3].append(NOISE_sig)
	"""
	domain = np.arange(len(sums)) 
	binsize = domain[1]-domain[0]
	plt.ylabel("Fraction of signal returned")
	plt.xticks(np.arange(1,len(sums)+1,step=1), ("6-order poly", "RCHPF", "Savitzky-Golay filter"), horizontalalignment='right', rotation=30)
	[plt.vlines(x=domain[i]+1, ymin=0, ymax=sums[i], linestyle = 'solid', linewidth = binsize*25, color = 'grey') for i in domain]
	[plt.text(x=domain[i]+1, y=sums[i]*(1+1/24),s="{0:0.2f} Percent".format(sums[i]*100)) for i in domain]
	plt.tight_layout()
	plt.savefig(savedir+"sums.pdf", dpi=600)
	plt.clf()
	"""
	#calculate dispersion error
	errors = numpy.asarray([numpy.std(i)/numpy.sqrt(len(i)) for i in signals])

	fig, ax = plt.subplots(1)
	fig.patch.set_visible(False)
	ax.axis('off')
	ax.axis('tight')
	ax.set_title("Dispersion Error")
	data = {"Sample Scan":[errors[0]], "RCHPF":[errors[1]], "6-order Poly":[errors[2]], "Savitzk Golay":[errors[3]]}
	df = pd.DataFrame(data = data)
	ax.table(cellText=df.values, colLabels=df.columns, loc='center')
	fig.tight_layout()

	plt.savefig(savedir+"dispersion_error.pdf", dpi=600)


	#Calculate gaussian histogram of autocorrelations
	plt.clf()
	plt.close('all')

	x = numpy.asarray(autocorrs_to_plot)
	means = means
	stdevs = stdevs
	for i in numpy.arange(len(x)):
		n, bins, patches = plt.hist(x[i], nbins, facecolor='b', alpha=0.5)
		
		text = r'$\mu={0:0.5f},\ \sigma={1:0.5f}$'.format(means[i],stdevs[i])
		plt.text(0.55,0.845, s=text, horizontalalignment='left', verticalalignment='baseline', transform=fig.transFigure)
		plt.xlabel('autocorrelation values')
		plt.ylabel('number of occurences')
		plt.title(autocorrs_titles[i])
		
		plt.savefig(savedir + "HIST" + autocorrs_filenames[i], dpi=600)
		plt.clf()
		
		
	#Run analysis with signal-less spectra
	"""
	BS_scan_res,_ = reciprocated_clone_hpf(sample_noise,3,**kwargs)
	BS_scan_bgless = BS_scan_res["filtereddata"].real
	SG_scan_bgless = SG(sample_noise)
	Poly_scan_bgless = poly_fit(sample_noise)
	
	sg_bgless_autocorrs = numpy.correlate(SG_scan_bgless-1., SG_scan_bgless-1., 'full')
	bs_bgless_autocorrs = numpy.correlate(BS_scan_bgless-1., BS_scan_bgless-1., 'full')
	poly_bgless_autocorrs = numpy.correlate(Poly_scan_bgless-1., Poly_scan_bgless-1., 'full')
	noise_autocorrs = numpy.correlate(sample_noise, sample_noise, 'full')
	"""
def controller(n, **kwargs):
	
	for i in numpy.arange(n):
		kwargs['iter']=i+1
		analysis(**kwargs)
		print("{0} of {1} runs complete".format(i+1, n))
	
	sig_atten = numpy.delete(sig_atten_array, 0, axis=1)
	sig_atten=[[eval(i) for i in x] for x in sig_atten]
	SG_sig_atten_avg = numpy.mean(numpy.asarray(sig_atten[0]))
	BS_sig_atten_avg = numpy.mean(numpy.asarray(sig_atten[1]))
	POLY_sig_atten_avg = numpy.mean(numpy.asarray(sig_atten[2]))
	NOISE = numpy.mean(numpy.asarray(sig_atten[3]))
	
	dispersion = numpy.delete(dispersion_array, 0, axis=1)
	dispersion = [[eval(i) for i in x] for x in dispersion]
	SG_dispersion_avg = numpy.mean(numpy.asarray(dispersion[0]))
	BS_dispersion_avg = numpy.mean(numpy.asarray(dispersion[1]))
	POLY_dispersion_avg = numpy.mean(numpy.asarray(dispersion[2]))
	NOISE_dispersion_avg = numpy.mean(numpy.asarray(dispersion[3]))
	
	SG_dispersion_error_avg = SG_dispersion_avg/nbins
	BS_dispersion_error_avg = BS_dispersion_avg/nbins
	POLY_dispersion_error_avg = POLY_dispersion_avg/nbins
	NOISE_dispersion_error_avg = NOISE_dispersion_avg/nbins
	
	fig, ax = plt.subplots(1)
	fig.patch.set_visible(False)
	ax.axis('off')
	ax.axis('tight')
	ax.set_title("Average Dispersion Error")
	data = {"RCHPF":[BS_dispersion_error_avg], "6-order Poly":[POLY_dispersion_error_avg], "Savitzk Golay":[SG_dispersion_error_avg], "Noise":[NOISE_dispersion_error_avg]}
	df = pd.DataFrame(data = data)
	ax.table(cellText=df.values, colLabels=df.columns, loc='center')
	fig.tight_layout()

	plt.savefig(savedir+"average_dispersion_error.pdf", dpi=600)
	plt.clf()
	
	#Plot signal attenuation
	plotter(numpy.asarray(sig_atten[1]), title="RCHPF signal attenuation", savedir = savedir+"RCHPF_sig_atten.pdf")
	plotter(numpy.asarray(sig_atten[0]), title="Savitzky Golay Signal Attenuation", savedir=savedir+"SG_sig_atten.pdf")
	plotter(numpy.asarray(sig_atten[2]), title="Polynomial fit Signal Attenuation", savedir=savedir+"poly_sig_atten.pdf")
	
	
	#cumulatively vertically add signals together 
	SG_coadded_signal = []
	BS_coadded_signal = []
	POLY_coadded_signal = []
	NOISE_coadded_signal = []
	del power_spec_array[0][0]
	del power_spec_array[1][0]
	del power_spec_array[2][0]
	del power_spec_array[3][0]
	SG_coadded_signal.append(power_spec_array[0][0])
	BS_coadded_signal.append(power_spec_array[1][0])
	POLY_coadded_signal.append(power_spec_array[2][0])
	NOISE_coadded_signal.append(power_spec_array[3][0])
	for i in numpy.arange(1, len(power_spec_array[0])):
		summed_array = [SG_coadded_signal[(i-1)][x]+power_spec_array[0][i][x] for x in numpy.arange(nbins)]
		SG_coadded_signal.append(summed_array)
	for i in numpy.arange(1, len(power_spec_array[1])):
		summed_array = [BS_coadded_signal[(i-1)][x]+power_spec_array[1][i][x] for x in numpy.arange(nbins)]
		BS_coadded_signal.append(summed_array)
	for i in numpy.arange(1, len(power_spec_array[2])):
		summed_array = [POLY_coadded_signal[(i-1)][x]+power_spec_array[2][i][x] for x in numpy.arange(nbins)]
		POLY_coadded_signal.append(summed_array)
	for i in numpy.arange(1, len(power_spec_array[3])):
		summed_array = [NOISE_coadded_signal[(i-1)][x]+power_spec_array[3][i][x] for x in numpy.arange(nbins)]
		NOISE_coadded_signal.append(summed_array)
	
	BS_signal_cumul_atten = []
	SG_signal_cumul_atten = []
	POLY_signal_cumul_atten = []
	NOISE_signal_cumul_atten = []
	
	
	for i in range(n):
		BS_sig = numpy.sum(BS_coadded_signal[i][128:133])
		SG_sig = numpy.sum(SG_coadded_signal[i][128:133])
		POLY_sig = numpy.sum(POLY_coadded_signal[i][128:133])
		NOISE_sig = numpy.sum(NOISE_coadded_signal[i][128:133])
		axion_sig = numpy.sum(axion_signal_array[i][128:133])
		BS_signal_cumul_atten.append(BS_sig/((i+1)*(axion_sig)))
		SG_signal_cumul_atten.append(SG_sig/((i+1)*(axion_sig)))
		POLY_signal_cumul_atten.append(POLY_sig/((i+1)*(axion_sig)))
		NOISE_signal_cumul_atten.append(NOISE_sig/((i+1)*(axion_sig)))
	plotter(numpy.asarray(BS_signal_cumul_atten), title="Coadded RCHPF signal atten", savedir = savedir + "coadded_RCHPF_signal_atten.pdf")
	plotter(numpy.asarray(SG_signal_cumul_atten), title="Coadded SG signal atten", savedir = savedir +  "coadded_SG_signal_atten.pdf")
	plotter(numpy.asarray(POLY_signal_cumul_atten), title="Coadded POLY signal atten", savedir = savedir +  "coadded_poly_signal_atten.pdf")
	plotter(numpy.asarray(NOISE_signal_cumul_atten), title="Coadded NOISE signal atten", savedir = savedir +  "coadded_noise_signal_atten.pdf")
		
	
	BS_autocorr_cumul = []
	SG_autocorr_cumul = []
	POLY_autocorr_cumul = []
	NOISE_autocorr_cumul = []
	BS_autocorr_cumul.append(autocorr_array[1][0])
	SG_autocorr_cumul.append(autocorr_array[0][0])
	POLY_autocorr_cumul.append(autocorr_array[2][0])
	NOISE_autocorr_cumul.append(autocorr_array[3][0])
	for i in numpy.arange(1, len(autocorr_array[1])):
		summed_array = [BS_autocorr_cumul[(i-1)][x]+autocorr_array[1][i][x] for x in numpy.arange(nbins-1)]
		BS_autocorr_cumul.append(summed_array)
	for i in numpy.arange(1, len(autocorr_array[0])):
		summed_array = [SG_autocorr_cumul[(i-1)][x]+autocorr_array[0][i][x] for x in numpy.arange(nbins-1)]
		SG_autocorr_cumul.append(summed_array)
	for i in numpy.arange(1, len(autocorr_array[2])):
		summed_array = [POLY_autocorr_cumul[(i-1)][x]+autocorr_array[2][i][x] for x in numpy.arange(nbins-1)]
		POLY_autocorr_cumul.append(summed_array)
	for i in numpy.arange(1, len(autocorr_array[3])):
		summed_array = [NOISE_autocorr_cumul[(i-1)][x]+autocorr_array[3][i][x] for x in numpy.arange(nbins-1)]
		NOISE_autocorr_cumul.append(summed_array)
	SG_cumul_std_devs = [numpy.std(i) for i in SG_autocorr_cumul]
	BS_cumul_std_devs = [numpy.std(i) for i in BS_autocorr_cumul]
	POLY_cumul_std_devs = [numpy.std(i) for i in POLY_autocorr_cumul]
	NOISE_cumul_std_devs = [numpy.std(i) for i in NOISE_autocorr_cumul]
	plotter(numpy.asarray(SG_cumul_std_devs), title="Coadded SG autocorr sigma", savedir = savedir + "Coadded_SG_autocorr_sigmas.pdf")
	plotter(numpy.asarray(BS_cumul_std_devs), title="Coadded BS autocorr sigma", savedir = savedir + "Coadded_RCHPF_autocorr_sigmas.pdf")
	plotter(numpy.asarray(POLY_cumul_std_devs), title="Coadded POLY autocorr sigma",savedir = savedir + "Coadded_poly_autocorr_sigmas.pdf")
	plotter(numpy.asarray(NOISE_cumul_std_devs), title="Coadded NOISE autocorr sigma",savedir = savedir + "Coadded_noise_autocorr_sigmas.pdf")

	
	
	
	
	
	
	
	String = "\n\nAverages of Signal attenuation calculation: \n\n Savitzky Golay: {0:0.5f} \n RCHPF: {1:0.5f} \n 6-order Poly fit: {2:0.5f} \n Original Noise: {3:0.5f}   (This should be close to unity)\n\n".format(SG_sig_atten_avg, BS_sig_atten_avg, POLY_sig_atten_avg, NOISE)
	String_disp = "\n\nAverages of Dispersion Calculation: \n\n Savitzky Golay: {0:0.5f} \n RCHPF: {1:0.5f} \n 6-order Poly fit: {2:0.5f} \n Original Noise: {3:0.5f}\n\n".format(SG_dispersion_avg, BS_dispersion_avg, POLY_dispersion_avg, NOISE_dispersion_avg)
	String_disp_error = "\n\nAverages of dispersion error: \n\n Savitzky Golay: {0:0.10f} \n RCHPF: {1:0.10f} \n 6-order Poly fit: {2:0.10f} \n Original Noise: {3:0.10f} \n\n".format(SG_dispersion_error_avg, BS_dispersion_error_avg, POLY_dispersion_error_avg, NOISE_dispersion_error_avg)
	print(String,String_disp, String_disp_error)







###############     Garbage dump     ###############

#Stuff I dont need but may need later so wont delete because I am lazy

"""
parser = argparse.ArgumentParser(description='A Background Subtraction (BS) script for removing large, low frequency structure')
parser.add_argument('-s', '--scanNumber', action='store', dest='scan_number', default = [1]*50, help='A Digitizer scan number or a set of data to be passed through the BS-er', nargs='*')
parser.add_argument('-n', '--numberOfClones', action='store', dest='NumberOfClonePairs', default = 3, help='The number of clones to be generated in the script')
parser.add_argument('-ch', '--channel', action='store', dest='ch', default=1, help='The channel used in the scan')
args = parser.parse_args()
scan_number = [float(i) for i in args.scan_number]
NumberOfClonePairs = int(args.NumberOfClonePairs)
ch = int(args.ch)
"""


"""
def poly_fit(signal):
        domain = numpy.arange(len(signal))
        polyfit = numpy.polyfit(domain, signal, 6)
        FilteredSignal = signal/(numpy.poly1d(polyfit)(domain))
        return FilteredSignal
    
def SG(signal):
    W = 101
    d=6
    to_remove = SGfilter(signal, W, d)
    FilteredSignal = signal/to_remove
    return FilteredSignal
"""



"""
	struc_length = 2*(numpy.absolute(argmax-argmin))
	sigmamult_n = (2*numpy.pi)/(struc_length)*5*10**2
	nnarg = lambda arr: (10*numpy.std(arr)/numpy.median(arr))**4
	frac_size = lambda arr: delta(arr)/arg_delta(arr)
	test_arr = lambda arr: (10**12)*frac_size(arr)*variation(arr)
	try_func = lambda arr: 2/(1 + (1/(arr**6)))
	test_sigma_func = lambda arr: 2**7*try_func(nnarg(arr))

"""
"""
	vratio = lambda arr: (numpy.max(arr)-numpy.min(arr))/(numpy.max(arr)+numpy.min(arr))
	wratio = lambda arr:  numpy.abs(numpy.argmax(arr) - numpy.argmin(arr))/len(arr)
	sigmamult = lambda arr: (vratio(arr)/wratio(arr))*5*10**2
	new_factor = lambda arr: numpy.std(arr)/numpy.mean(arr)
	sigma_func = lambda arr, n: numpy.ceil((2*n+1)*sigmamult(arr))
	
	
	test_arg = lambda arr: ((numpy.max(arr)-numpy.min(arr))/numpy.std(arr))
	test_sigma_func = lambda arr: numpy.abs(numpy.tanh(test_arg(arr)))**0.85*(numpy.max(arr)-numpy.min(arr)) + (1-numpy.abs(numpy.tanh(test_arg(arr))))**0.85*len(arr)
"""

"""
	
	func = lambda a,b,c,arr: a*numpy.sin(b*(arr-c))/(b*(arr-c))
	p_guess = [numpy.std(ORIG_corr), numpy.pi/arg_delta(ORIG_corr), len(ORIG_corr)/2]
	popt, pcov = curve_fit(func, numpy.arange(len(ORIG_corr)), ORIG_corr, p0=p_guess)
	ORIG_corr_func = func(popt[0], popt[1], popt[2], numpy.arange(len(ORIG_corr)))
	width_factor = numpy.mean(arr)**(-2)*delta(ORIG_corr_func)
"""