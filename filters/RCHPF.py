import numpy
import astropy.constants as C
import astropy.units as u
from math import gamma
from scipy.signal import savgol_filter as SGfilter
#import argparse
import pyopencl as cl
import reikna.cluda as cluda
from reikna.fft import FFT as fftcl
import os
import sys
sys.path.append("..")
from toolbox.DFT import DFT, IDFT
import time

os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1"
####################################
savedir = "C:/Users/drums/Documents/ADMX/Papers/Own Papers/Plots and Graphics/Bsub/MyBSub/"
####################################


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


# data and axion params 
nsamp = 1
c = C.c.to(u.km/u.s).value # in km/s
h = C.h.to(u.eV*u.s).value # in eV*s
k = C.k_B.value #J/K
ma = 1.e-6
rest_energy = ma
nbins = 2**8 # Length of scan
df = 100.0 #Bin Size
start = rest_energy/h - nbins*df/2
stop = start + nbins*df
bins = [start + df*i for i in range(nbins)]


	
def reciprocated_clone_hpf(data, npairs, testing=False, scan_number='(Unknown)', **meta):
	"""
	Parameters:
				data: can be array of data or a digitizer number
				npairs: number of clones to make when generating large fourier domain
				(removed)ch: --OBSOLETE-- specify channel data was taken on.
				testing: Used for testing... Plots major calculations in the script of meta analysis if True
				scan_number: used for testing. Specifies scan number in title when saving plots
	"""
	
	avail_types = [numpy.ndarray, numpy.arange(1), 1, 1.]
	for i in avail_types:
		if isinstance(data, i):
			break
		else:
			return "Error: invalid data type"
	errors = {'maxed_filter_size':False}
	#if (numpy.mean(data)-data[0])>numpy.std(data): # Some spectra have odd behavior at beginning of scan, i.e. a single downward spike at the beginning position. I just set a default value
	data[0]=data[1]
	while len(data)<256: #For some reason, some datasets dont have 2**8 elements.
		data = numpy.append(data[0], data)
	
	reflecting_time_start = time.time()
	ORIG = data
	sz = len(data)
	n = npairs # Specifies number of concatenated datasets
	REFL = data[::-1] # Generate a reflected copy of data
	DOMAIN = numpy.arange(len(data)) # Generate domain for edge-matching function
	reflecting_time_stop = time.time()
	
	#Generate statistics about the data set and define variables

	calculating_tophat_size_start = time.time()
	maxdata = numpy.max(data)
	mindata = numpy.min(data)
	argmax = numpy.argmax(data)
	argmin = numpy.argmin(data)
	"""
	vratio = (maxdata-mindata)/(maxdata+mindata)
	wratio = numpy.absolute(argmax - argmin)/sz
	sigmamult = (vratio/wratio)*5*10**2
	sigmao = numpy.ceil((2*n+1)*sigmamult)
	"""
	
	struc_length = 2*(numpy.absolute(argmax-argmin))
	sigmamult = (2*numpy.pi)/(struc_length)*5*10**2
	sigma = numpy.ceil(sigmamult*(2*n+1))
	if sigma>200:
		sigma=200
		print("Ceiling of reciprocated clone hpf reached")
		errors['maxed_filter_size'] = True
		
	calculating_tophat_size_stop = time.time()
	
	#reflects the dataset about the center and operates on it to sequence datasets without discontinuities 
	
	reciprocating_array_start = time.time()
	a = data[0]
	b = data[-1]
	RECIP = ((b**2 - a**2)*(numpy.cos(DOMAIN*numpy.pi/(2*numpy.size(REFL)))**2) + a**2)/REFL
	reciprocating_array_stop = time.time()
	
	#Construct large array consisting of reflected and nonreflected data_array
	generating_large_array_start = time.time()
	SigToPlot = [*RECIP, *ORIG]*(n+1)
	szA = len(SigToPlot) # Find length of large array
	generating_large_array_stop = time.time()
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
	fftSigToPlot = DFT(SigToPlot) 
	reducedfftSigToPlot = tophat*fftSigToPlot # Isolate lower frequencies
	fft_and_highpassfilter_stop = time.time()
	
	
	# Inverse Fourier Transform
	ifft_start = time.time()
	BSfftSigToPlot = IDFT(reducedfftSigToPlot)
	ifft_stop = time.time()
	
	#pick out the original signal
	
	picking_og_signal_start = time.time()
	pickORIG = [0]*len(data)
	if n%2==1: #calculate which part of the array to pick out. if n is odd, pick out the original scan to the left of the center of the array. if n is even, pick out the middle
		l=n
	elif n%2==0:
		l=n+1
		
	for i in range(len(data)):
		pickORIG[i] = BSfftSigToPlot[(l*len(data)+i)]
	picking_og_signal_stop = time.time()
	
	dividing_structure_start = time.time()
	ORIG = numpy.array(data) #Convert back into an array
	pickORIG =  numpy.array(pickORIG)
	filtereddata = ORIG/pickORIG # Divide out low freq. structure from data.
	dividing_structure_stop = time.time()
	
	if testing == True:
		from toolbox.plot_dataset import plotter
		sdir = "C:/Users/drums/Documents/Coding Software/Python/Scripts/New-Analysis-Scheme/oo_analysis/figures/"+ '('+scan_number+')'
		format = '.pdf'
		print("Plotting data")
		plotter(data,savedir = sdir+'data'+format); 
		print("Plotting reflected data");
		plotter(REFL,savedir = sdir+'refl_data'+format); 
		print("Plotting reciprocated array");
		plotter(RECIP,savedir = sdir+'recip_datas'+format); 
		print("Plotting large appended array")
		plotter(SigToPlot,savedir = sdir+'large_append_array'+format); 
		print("Plotting fft'ed large appended array");
		plotter(fftSigToPlot,savedir = sdir+'fft-ed large appended_array'+format, **{'title': "With sigma = {0}".format(sigma)}); 
		print("Plotting inverse fft'ed subtracted fft")
		plotter(BSfftSigToPlot,savedir = sdir+'inverse fft-ed append array'+format);
		print("Plotting recovered data from inverse fft'ed array")
		plotter(pickORIG,savedir = sdir+'recovered data from ifft-ed array'+format)
		print("Plotting filtered data")
		plotter(filtereddata,savedir = sdir+'filtered data'+format)
		print("Plotting filter function")
		plotter(tophat, savedir = sdir + "filter_function"+format)
	
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
	
	return {"filtereddata":filtereddata, "sigma": sigma, "number of clones": npairs, "meta": meta}, errors

def sixorderpoly(signal):
        domain = numpy.arange(len(signal))
        polyfit = numpy.polyfit(domain, signal, 6)
        FilteredSignal = signal/(numpy.poly1d(polyfit)(domain))
        return FilteredSignal
    
def SG_filter(signal):
    W = 101
    d=6
    to_remove = SGfilter(signal, W, d)
    FilteredSignal = signal/to_remove
    return FilteredSignal

def ADMX_signal_model(freq,exp_pow,poly_pow,temp):
    if freq<=rest_energy/h:
        return 0.0
    else:
        pass
    sig = temp*rest_energy/h
    nuosig = (freq-rest_energy/h)/sig
    gammafunc = gamma((1.0+poly_pow)/exp_pow)
    Cnum = exp_pow
    Cden = ((1.0/(sig))**exp_pow)**(-(1.0+poly_pow)/exp_pow)*(1.0/sig)**(poly_pow)*gammafunc
    C = Cnum/Cden
    dist = C*(nuosig)**poly_pow*numpy.exp(-(nuosig)**exp_pow)
    return dist
signal_strength = 0.15/3 # Level to be recoverable with ~20 scans at SNR=3
params = [1.4,0.37,4.7e-7]
#plotter(numpy.array([signal_strength*df*ADMX_signal_model(x, params[0], params[1], params[2]) for x in bins])[128:134])
	
def generate_signal(with_signal=True):
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
	if with_signal:
		Combined_Signal = numpy.array((noise+signals+1.0)*background)[0]
	else:
		Combined_Signal = numpy.array((noise+1.0)*background)[0]
	noises = noise+signals
	toreturn = {"Combined signal":Combined_Signal, "noise":noises[0], "backgrounds":background[0], "axion signal":single_signal}
	return toreturn
""" 
def plotter(signal, title=None, xlabel=None, ylabel=None, Min = None, Max = None, filename="None"):
	from matplotlib import pyplot as plt
	#import pandas as pd
	plt.style.use('seaborn-pastel')
	#plt.clf() #first clear current figures, as safe guard

	if Min==None:
		Min=numpy.min(signal)
	if Max==None:
		Max=numpy.max(signal)
		
	f, ax = plt.subplots(1, figsize=(6,6))
	length = len(signal)
	
	ax.set_title(title)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_ylim(Min, Max)
	savefile = savedir + filename
    
	ax.plot(numpy.arange(length), signal.real, linestyle='-')
	plt.tight_layout()

	if filename=="None":
		plt.show()
	else:
		plt.savefig(savefile, dpi=600)
	#plt.clf() #clear figures again
"""	


	
	

	
	
	
	
	
	
	
	
	
