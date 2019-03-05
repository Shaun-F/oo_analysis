import sys; sys.path.append('..')
import numpy as np
import toolbox.lorentzian
from toolbox.lorentzian import lorentzian_value

def axion_power(axion_dataset, frequency, form_factor=0.4, axion_coupling=0.36, halo_density=0.45, nbins=2**8, res=100):
	#Halo density in Gev/(c^2 cm^3)
	#frequency in Hz
		
	scan = axion_dataset
	if scan=="error" or scan=='None':
		return "error: Fault in axion dataset fetching"
	
	#extract relevent experimental parameters
	q = scan.attrs["Q"]
	B = scan.attrs["b_field"]
	vol = 220 
	res_freq = (scan.attrs["mode_frequency"])*10**6 # in Hz
	l_value = lorentzian_value(q, res_freq, frequency)
	
	total_power = (2.09*10**(-22))*(vol/220)*((B/7.6)**2)*form_factor*((axion_coupling/0.36)**2)*(halo_density/0.45)*(res_freq/750)*(q/70000)*l_value
	
	return total_power
	
	
	