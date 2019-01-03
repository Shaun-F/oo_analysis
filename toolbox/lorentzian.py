import numpy as np

def lorentzian(Q, resonant_frequency, nbins, resolution):
	freqs = np.arange(start = (resonant_frequency-(resolution*nbins/2)), stop = (resonant_frequency + (resolution*nbins/2)), step = resolution)

	l = lambda FWHM, center, f: (1/np.pi)*(0.5*Q)/((f-center)**2 + (0.5*Q)**2)

	return l(Q, resonant_frequency, freqs)