import numpy as np 
import astropy.constants as C
import astropy.units as u
from math import gamma
from scipy.signal import savgol_filter as SGfilter
#import argparse
import pyopencl as cl
import reikna.cluda as cluda
from reikna.fft import FFT as fftcl
import os

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

def DFT(inputsignal):
    signal = np.array(inputsignal, dtype = np.complex128)
    device = cl.get_platforms()[1].get_devices(device_type=cl.device_type.GPU)
    ctx = cl.Context(device)
    queue = cl.CommandQueue(ctx)
    
    api = cluda.ocl_api()
    thr = api.Thread(queue)
    signal_dev = thr.to_device(signal)
    fft_res_dev = thr.array(signal.shape, dtype = np.complex128)
    
    FFT = fftcl(signal_dev).compile(thr)
    FFT(fft_res_dev, signal_dev)
    
    res = fft_res_dev.get()
    return res

def IDFT(inputsignal):
    
    signal = np.array(inputsignal, dtype = np.complex128)
    device = cl.get_platforms()[1].get_devices(device_type=cl.device_type.GPU)
    ctx = cl.Context(device)
    queue = cl.CommandQueue(ctx)
    
    api = cluda.ocl_api()
    thr = api.Thread(queue)
    
    signal_dev = thr.to_device(signal)
    ifft_res_dev = thr.array(signal.shape, dtype = np.complex128)
    
    IFFT = fftcl(ifft_res_dev).compile(thr)
    IFFT(ifft_res_dev, signal_dev, inverse = True)
    
    res = ifft_res_dev.get()
    return res
	
def reciprocated_clone_hpf(data, npairs):
	"""
	Parameters:
				data: can be array of data or a digitizer number
				npairs: number of clones to make when generating large fourier domain
				ch: --OBSOLETE-- specify channel data was taken on.
	"""
	avail_types = [type(np.arange(1)), type(1), type(1.)]
	if type(data) in avail_types:
		pass
	elif type(data) not in avail_types:
		return "Error: invalid data type"
	
	ORIG = data
	sz = len(data)
	n = npairs # Specifies number of concatenated datasets
	REFL = data[::-1] # Generate a reflected copy of data
	DOMAIN = np.arange(len(data)) # Generate domain for edge-matching function

	
	#Generate statistics about the data set and define variables

	maxdata = np.max(data)
	mindata = np.min(data)
	argmax = np.argmax(data)
	argmin = np.argmin(data)
	vratio = (maxdata-mindata)/(maxdata+mindata)
	wratio = np.absolute(argmax - argmin)/sz
	sigmamult = (vratio/wratio)*5.0+3.0
	sigma = np.ceil((2*n+1)*sigmamult)


	#reflects the dataset about the center and operates on it to sequence datasets without discontinuities 
	a = data[0]
	b = data[-1]
	RECIP = ((b**2 - a**2)*(np.cos(DOMAIN*np.pi/(2*np.size(REFL)))**2) + a**2)/REFL

	#Construct large array consisting of reflected and nonreflected data_array
	SigToPlot = [*RECIP, *ORIG]*(n+1)
	szA = len(SigToPlot) # Find length of large array
	tophat = [0]*szA # Generate tophat function to pick out lower frequency structure
	for i in range(int(sigma)):
		tophat[i]=1
		tophat[szA-i-1] = 1
	
	
	# Fourier transform large array
	fftSigToPlot = DFT(SigToPlot) 
	reducedfftSigToPlot = tophat*fftSigToPlot # Isolate lower frequencies
	
	
	
	# Inverse Fourier Transform
	BSfftSigToPlot = IDFT(reducedfftSigToPlot)
	
	#pick out the original signal
	pickORIG = [0]*len(data)
	for i in range(len(data)):
		pickORIG[i] = BSfftSigToPlot[(3*len(data)+i)]
		
	ORIG = np.array(data) #Convert back into an array
	pickORIG =  np.array(pickORIG)
	filtereddata = ORIG/pickORIG # Divide out low freq. structure from data.
	"""
	teller=0
	for i, val in enumerate(filtereddata):
		if val.imag<10**(-3):  #NOTE: this restriction on the imaginary part may be too lenient. 
			teller += 0
		else:
			teller += 1
	if teller==0:
		filtereddata = filtereddata.real
	elif teller>0:
		return (print("\n\n\nERROR: Returned signal of RCHPF should be real\n\n\n"), print(filtereddata), print(sigma))
	"""
	return {"filtereddata":filtereddata, "sigma": sigma, "number of clones": npairs}

def sixorderpoly(signal):
        domain = np.arange(len(signal))
        polyfit = np.polyfit(domain, signal, 6)
        FilteredSignal = signal/(np.poly1d(polyfit)(domain))
        return FilteredSignal
    
def SG_filter(signal):
    W = 101
    d=6
    to_remove = SGfilter(signal, W, d)
    FilteredSignal = signal/to_remove
    return FilteredSignal
    
def plotter(signal, title=None, xlabel=None, ylabel=None, Min = None, Max = None, filename="None"):
	from matplotlib import pyplot as plt
	#import pandas as pd
	plt.style.use('seaborn-pastel')
	#plt.clf() #first clear current figures, as safe guard

	if Min==None:
		Min=np.min(signal)
	if Max==None:
		Max=np.max(signal)
		
	f, ax = plt.subplots(1, figsize=(6,6))
	length = len(signal)
	
	ax.set_title(title)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_ylim(Min, Max)
	savefile = savedir + filename
    
	ax.plot(np.arange(length), signal.real, linestyle='-')
	plt.tight_layout()

	if filename=="None":
		plt.show()
	else:
		plt.savefig(savefile, dpi=600)
	#plt.clf() #clear figures again
	
def ADMX_signal_model(freq,exp_pow,poly_pow,temp):
    if freq<=rest_energy/h:
        return 0.0
    else:
        pass
    sig = temp*rest_energy/h
    nuosig = (freq-rest_energy/h)/sig
    gamma = gamma((1.0+poly_pow)/exp_pow)
    Cnum = exp_pow
    Cden = ((1.0/(sig))**exp_pow)**(-(1.0+poly_pow)/exp_pow)*(1.0/sig)**(poly_pow)*gamma
    C = Cnum/Cden
    dist = C*(nuosig)**poly_pow*np.exp(-(nuosig)**exp_pow)
    return dist
signal_strength = 0.15/3 # Level to be recoverable with ~20 scans at SNR=3
params = [1.4,0.37,4.7e-7]
#plotter(np.array([signal_strength*df*ADMX_signal_model(x, params[0], params[1], params[2]) for x in bins])[128:134])
	
def generate_signal(with_signal=True):
	signal_strength = 0.05/3 # Level to be recoverable with ~20 scans at SNR=3
	params = [1.4,0.37,4.7e-7]
	order = 3
	#random polynomials
	coeffs = (2*np.random.random((order,))-1.0)
	rand_poly = lambda x: sum([coeff*x**i for i,coeff in enumerate(coeffs)])
	#random hot rod
	Q = np.random.uniform(low=10000, high=60000)
	coeff = np.random.uniform(low=0.0, high=1)
	rand_hot_rod = lambda x: coeff/(1.0+(x-rest_energy/h)**2/(rest_energy/(h*Q))**2)
	
	polys=[]
	for i in range(nsamp):
		polys.append([0.1*rand_poly(1.0*x/nbins) for x in range(nbins)])
		
	ress = []
	for i in range(nsamp):
		ress.append([rand_hot_rod(x) for x in bins])
	background = np.ones((nsamp, nbins)) + np.array(polys) + np.array(ress)
	
	single_signal = np.array([signal_strength*df*ADMX_signal_model(x, params[0], params[1], params[2]) for x in bins])
	signals = np.array([single_signal for i in range(nsamp)])
	#Create noise test signals
	noise = 0.01*np.random.rand(nsamp, nbins)
	#combine everything
	if with_signal:
		Combined_Signal = np.array((noise+signals+1.0)*background)[0]
	else:
		Combined_Signal = np.array((noise+1.0)*background)[0]
	noises = noise+signals
	toreturn = {"Combined signal":Combined_Signal, "noise":noises[0], "backgrounds":background[0], "axion signal":single_signal}
	return toreturn
 


	
	

	
	
	
	
	
	
	
	
	
