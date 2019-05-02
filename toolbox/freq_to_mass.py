import astropy.constants as C
import astropy.units as U

def freq_to_mass(freq):
	"""
	function takes in a frequency in units of Hz and outputs the corresponding mass in units of eV
	"""
	return (freq*(1/U.s)*C.h).to(U.eV).value