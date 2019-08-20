import numpy as np

def lorentzian(Q, resonant_frequency, starting_freq, stoping_freq, resolution):
	"""
	Calculate lorentzian function with starting frequency determined by starting_freq, ending frequency determined by stoping_freq, binwidth determined by resolution,
	central frequency determined by resonant_frequency, and width determined by Q
	"""
	if resonant_frequency<10**3:
		resonant_frequency*=10**6 #force to Hz
	if resolution<1:
		resolution*=10**6 #Force to Hz
	if starting_freq<1000:
		starting_freq*=10**6
	if stoping_freq<1000:
		stoping_freq*=10**6
	
	freqs = np.arange(start = starting_freq, stop = stoping_freq, step = resolution)
	
	
	l = lambda Q, center, f: 1/(1 + 4*((f-center)**2)*(Q**2/center**2))    #FWHM = resonant_frequency/Quality_factor
	
	return l(Q, resonant_frequency, freqs)
	
def lorentzian_value(Q, resonant_frequency, freq):
	
	l = lambda Q, center, f: 1/(1 + 4*((f-center)**2)*(Q**2/center**2))

	return l(Q, resonant_frequency, freq)

