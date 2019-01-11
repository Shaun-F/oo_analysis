import numpy as np

def lorentzian(Q, resonant_frequency, nbins, resolution):
	freqs = np.arange(start = (resonant_frequency-(resolution*nbins/2)), stop = (resonant_frequency + (resolution*nbins/2)), step = resolution)
	
	#Force length 
	if len(freqs)<nbins:
		while len(freqs)<nbins:
			freqs = np.append(freqs, freqs[-1])
	elif len(freqs)>nbins:
		while len(freqs)>nbins:
			freqs = np.delete(freqs, -1)
			
	l = lambda FWHM, center, f: (1/np.pi)*(0.5*Q)/((f-center)**2 + (0.5*Q)**2)

	return l(Q, resonant_frequency, freqs)
	
def lorentzian_value(Q, resonant_frequency, freq):
	
	l = lambda f: (1/np.pi)*(0.5*Q)/((f-resonant_frequency)**2 + (0.5*Q)**2)
	return l(freq)