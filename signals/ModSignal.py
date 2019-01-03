import numpy as np
import scipy as sp
from astropy import constants as Cnsts
from datetime import datetime as dt
import pytz
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.pylab as plt
import pandas as pd
import argparse

"""
Could make this script easier to run with argument parsing
"""
parser = argparse.ArgumentParser()
parser.add_argument('-s','--signal', action='store', default = 'SHM', dest = 'signal_model')
args = parser.parse_args()

shape_model = args.signal_model



def timestamptranslator(timestamp):
	"""
	Description:Function takes in a time stamp either in the form 'dd-mm-yyyy  hh:mm:ss AM' or 'yyyy-mm-dd hh:mm:ss-ss' and returns the total hours elapsed since June 17th, 1998 midnight GMT (Input timezone is PST).
	Parameters: timestamp
	Output: Number of hours since Midnight june 17th, 1998 GMT from input timestamp

	"""

	if timestamp[2]== '-':
		"""This IF statement tests whether the date string is in the form 'dd-mm-yyyy  hh:mm:ss AM'. This is the case for a digitizer timestamp"""
		d1 = dt.strptime(timestamp, "%d-%m-%Y  %I:%M:%S %p")
		d1 = dt(d1.year,d1.month,d1.day,d1.hour,d1.minute,d1.second,0, tzinfo=pytz.timezone('America/Vancouver'))
	if timestamp[4]== '-':
		"""This IF statement tests whether the date string is in the form 'yyyy-mm-dd hh:mm:ss-ss'. This is the standard ISO 8601 format and the format Python uses for the datetime module"""
		d1 = dt.strptime(timestamp, "%Y-%m-%d %H:%M:%S-%f")
		d1 = dt(d1.year,d1.month,d1.day,d1.hour,d1.minute,d1.second,0, tzinfo=pytz.timezone('America/Vancouver'))
	tz_pst = pytz.timezone('America/Vancouver') #Pacific standard time
	tz_grnwch = pytz.timezone('Etc/Greenwich') #greenwich mean time

	d2 = dt(1949,12,31,23,58,56,3795, tzinfo = tz_grnwch) #B1950 epoch in UTC, the standard UNIX reference date
	return ((d1-d2).total_seconds()/(60*60)) #Number of hours since B1950

def Equ_to_Gal(numberofhours):
	"""
	Description: Function returns the velocity vector of the lab in galactic cylindrical coordinates using a matrix equation found in Spherical Astronomy by Robin Green
	Parameters: Number of hours since B1950
	Output:Velocity vector in Galactic cylindrical Coordinates
	"""

	phi = 123 # angle between vr and vt at the start of B1950 epoch
	v=0.4638*np.cos(47.66) #Velocity of experiment in equatorial frame. 0.4638 is rotation speed of earth, 47.66 is latitude of experiment

	t = numberofhours*3600 # sidereal day in seconds after midnight
	T = 86164.09056 # a sidereal day in seconds
	vx = v*np.cos(2*np.pi*t/(T+phi))
	vy = v*np.sin(2*np.pi*t/(T+phi))
	vz = 0

	vr = -0.0669887*vx + -0.8727558*vy + -0.4835389*vz
	vt = 0.4927285*vx + -0.4503470*vy + 0.7445846*vz
	vz  = -0.8676008*vx + -0.1883746*vy + 0.4601998*vz
	#Takes in the velocity vector of the experiment in equatorial cartesian basis and transforms to the galactic cylindrical basis

	return {"vr":vr, "vt":vt, "vz":vz}

def earth_orbit_velocity(timestamp):
	"""
	Description: Function returns the orbital velocity in the cylindrical galactic center basis vr, vt, vz
	Parameters:Timestamp
	Output: Dictionary of velocity components 'vr', 'vt', 'vz'
	"""
	time = timestamptranslator(timestamp)
	t = time*86164.09056/86400
	T = 3.1558*10**7
	phi = 86.14
	v = 29.785
	v_plane = {"vx":v*np.cos(2*np.pi*t/(T+phi)), "vy":v*np.sin(2*np.pi*t/(T+phi)), "vz":0}
	v_gal = {"vr":Equ_to_Gal(time)["vr"]*v_plane["vx"], "vt":Equ_to_Gal(time)["vt"]*v_plane["vy"]}

	return v_gal

def cavity_rotation_velocity(timestamp):
	"""
	Description: Function returns the cavitys orbital velocity in cylindrical galactic center basis vr, vt, vz
	Parameters: timestamp
	Output: Dictionary to velocity components 'vr', 'vt', 'vz'
	"""
	time = timestamptranslator(timestamp)
	t = time*3600
	T = 86164.09056
	phi = -122.3
	v = 0.4638*np.cos(47.66)
	v_spin = {'vx': v*np.cos(2*np.pi*t/(T+phi)), 'vy':v*np.cos(2*np.pi*t/(T+phi)), 'vz':0}
	v_gal = {'vy': Equ_to_Gal(time)['vr']*v_spin['vx'], 'vt': Equ_to_Gal(time)['vt']*v_spin['vy']}

	return v_gal

def modulation(modulation_type, timestamp):
	date = timestamp
	hourtime = timestamptranslator(date)

	vel_orbit = (earth_orbit_velocity(date))
	vt_orbit = vel_orbit["vt"]

	vel_rotation = cavity_rotation_velocity(date)
	vt_rotation = vel_rotation["vt"]

	if modulation_type == None or modulation_type =="None":
		modulation_type = "solar"

	if modulation_type == "solar":
		return 0
	elif modulation_type=="earth orbit":
		return vt_orbit
	elif modulation_type=="earth rotation":
		return vt_orbit + vt_rotation
	else:
		return "Error: modulation type not recognized"

def axionDM_w_baryons(mod_vel, v_sol, mass, freq):

	"""
	Description: Function gives probability density for axion signal of Lentz et al. 2017
	Parameters: mod_vel (km/s) is the variation of the experiments velocity from the mean solar velocity w.r.t a right-handed galactic coordinate system centered on the galactic center. v_sol (km/s) is the mean solar velocity. mass (eV) is the mass of the axion. freq (MHz) is a desired sample point.
	Output: Probability to find an axion with given input frequency.
	"""
	# Isolate the tangential velocity
	if type(mod_vel)==type(0.0):
		v_m = mod_vel
	elif mod_vel.tangential:
		v_m = mod_vel.tangential
	elif mod_vel[1]:
		v_m = mod_vel[1]
	else:
		return "Error: Modulation velocity format not recognized"


	f = freq*10**6 #Hz
	c = Cnsts.c.value*10**(-3) #Km/s
	RME = mass #eV
	h = Cnsts.h.value*6.242*10**18 #eV*s
	f_o = RME/h #rest mass frequency in 1/s

	if f_o>f: #If frequency is less than rest mass frequency, probability is zero.
		return 0.0

	#Define various parameters and constants for the distribution function
	alpha = 0.3
	beta = 1.3
	T = 4.7*10**(-7)
	gamma = sp.special.gamma((1+alpha)/beta)

	Boost = 1-(mod_vel)*(v_sol)/(v_sol**2) #Modulating parameter to account for diurnal variations



	N_num = beta
	N_den = ((Boost/(T*f_o))**beta)**(-(1.0+alpha)/beta)*(Boost/(T*f_o))**(alpha)*gamma
	nrmcnst = N_num/N_den

	#The distribution function
	A = nrmcnst*(((f-f_o)*Boost)/(f_o*T))**(alpha)*np.exp(-(((f-f_o)*Boost)/(f_o*T))**beta)

	return A

def IsThmSph(freq,mass,sigma_v,v_lab):
	"""Description: Function takes in parameters and outputs probability to find axion at input frequency
	Parameters: frequency (MHz) to find corresponding probability. Mass of axion in ev. velocity dispersion of axions in Kilometers per second. velocity of laboratory in galactic frame in Kilometers per second.
	Output: Probability to find axion at input frequency"""

	c = Cnsts.c.value*10**(-3) #Speed of light in kilometers per second
	h = Cnsts.h.value*6.242*10**18 #Plancks constant in units eV seconds
	RME = mass #Rest Mass energy of Axion in eV
	m = RME/(c**2) #Mass of Axion in eV/c^2
	E = freq*h*10**(6) #Energy of equivalent photon at input frequency in Hz
	KE = E-RME #Kinetic energy of axion in eV
	rmfreq = (RME/h)*10**(-6) # Rest mass frequency of axion
	v = c*np.sqrt((2*KE)/RME) #velocity of axion
	beta = 1/(2*(sigma_v**2)) #this is the beta from turner 1990


	if rmfreq<freq:
		X = (2*h*(c**2)*np.sqrt(beta/np.pi))/(RME*v_lab)
		Y = np.exp(-beta*(v_lab**2+v**2))
		Z = np.sinh(2*beta*v_lab*v)
		dist = X*Y*Z
	else: #If input frequency is less than the rest mass frequency, probability is zero (nonphysical)
		dist = 0

	return dist

def SHM(freq, mass, modulation_vel):
	sigma_v = 160
	v_lab = modulation_vel + sigma_v*np.sqrt(2)

	return IsThmSph(freq, mass, sigma_v, v_lab)

def pCDM_maxwell_like_form(mass,freq,v_lab, alpha,beta,T):
	"""
	Description:Function gives the frequency power spectra for halo axions as modeled by the N-Body simulation using pCDM, parameterized by three values alpha, beta, T
	Parameters:mass of axion (eV), sampling frequency (Hz), velocity of laboratory (km/s), parameter values of N-Body simulation
	Output: Power spectra of axion
	"""
	f = freq*10**6 #sampling frequency
	c = Cnsts.c.value*10**(-3) #Speed of light
	RME = mass #Rest mass energy of axion
	h = Cnsts.h.value*6.242*10**18 # plancks constant in ev*s
	fo = RME/h # Rest mass frequency of axion
	vsol = 240 # average solar speed in km/s

	gamma = sp.special.gamma((1.0+alpha)/beta)
	Cnum = beta
	Boost = 1-v_lab*vsol/(vsol**2)
	Cden = ((Boost/(T*fo))**beta)**(-(1.0+alpha)/beta)*(Boost/(T*fo))**(alpha)*gamma
	Cnst = Cnum/Cden

	if fo<f:
		A = Cnst*(h*(f-fo)/(mass*T))**(alpha)*np.exp(-(h*(f-fo)/(mass*T))**beta)
	else:
		A = 0

	return A

def pCDM_only(mass,freq, v_lab):
	"""
	Description: Function returns power spectra for pressureless cold dark matter without baryons
	Parameters:mass of axion (eV), Sampling frequency (Hz)
	Output: Power spectra at sampling frequency
	"""
	alpha = 0.25
	beta = 1.3
	T = 2.0*10**(-7)

	return pCDM_maxwell_like_form(mass,freq, v_lab, alpha,beta,T)

def pCDM_w_baryons(mass,freq, v_lab):
	"""
	Description: Function returns power spectra for pressureless cold dark matter, taking into account baryons
	Parameters: Mass of axion (eV), Sampling frequency (Hz)
	Output: Power spectra at sampling frequency
	"""
	alpha = 0.36
	beta = 1.39
	T = 4.7*10**(-7)

	return pCDM_maxwell_like_form(mass,freq,v_lab,alpha,beta,T)


def modulatedsignal(modulation_type, timestamp, shape_model, axion_mass, resolution, wantTseries, startfreq, alpha, beta, T):
	"""
	Description:Function bins a given velocity-modulated axion shape by integration.

	Parameters: Modulation types include 'solar', 'earth orbit', 'earth rotation'. Timestamp takes either the form "dd-mm-yyyy  hh:mm:ss AM", which is used by the Digitizer, or the form "yyyy-mm-dd hh:mm:ss-ss", which is the standard ISO 8601 format. Current shape models include "SHM", "axionDM_w_baryons", "pCDM_only", "pCDM_w_baryons", "pCDM_maxwell_like_form", "axionDM_only". Axion mass is in eV.

	Optional Parameters: Resolution give the size of the bins in Hz. wantTseries specifies if the output should be in frequency space or as a timeseries (y/n). startfreq allows you to specify the position of the signal shape within the array (MHz). alpha, beta, T are fit parameters for the N-Body model.

	Output: velocity-modulated signal shape over frequency or time.
	"""


	if type(resolution) != float:
		res = 100
	if type(wantTseries) != str:
		wantTseries='n'
	if type(alpha) != float:
		alpha = 0
	if type(beta) != float:
		beta = 0
	if type(T) != float:
		T = 0

	#specify parameters to be used in various signal scripts
	modtype = modulation_type #create pointer with shorter character length
	vsol = 232 # (km/s) mean solar speed around galaxy
	vel_mod = float(modulation(modtype, timestamp)) #variation of experimental speed around mean solar speed
	h = Cnsts.h.value*6.242*10**18 # plancks constant in eV*s
	m = axion_mass # axion mass in eV
	bin_width = float(resolution) #size of bins in hZ
	RMF = float(m/h) # Rest mass frequency in hZ
	cutoffarea = 0.9


	#set default value for starting frequency
	if type(startfreq) != float and type(startfreq) != int:
		startfreq = RMF

	binsize = []
	binned_signal = []
	sumbins = 0

	startingfreqs = []
	info = []
	i=0


	#modulate specified axion shape based on input parameters
	if shape_model == "axionDM_w_baryons":
		while sumbins<cutoffarea:
			scanfreq = float(startfreq + (i)*(bin_width/3))*10**(-6)
			binsize.append(axionDM_w_baryons(vel_mod, vsol, m, scanfreq)*(bin_width/3))
			if ((i+1)/3) == np.floor((i+1)/3) and i>0:
				binned_signal.append(binsize[int(i-2)] + binsize[int(i-1)] + binsize[int(i)])
				sumbins = sumbins + binned_signal[int(i/3)]
				startingfreqs.append(scanfreq-bin_width)
			i+=1
			if len(binned_signal)>500:
				break

	elif shape_model == "SHM":
		while sumbins<cutoffarea:
			scanfreq = float(startfreq+(i)*(bin_width/3))*10**(-6)
			binsize.append(SHM(scanfreq, m, vel_mod)*(bin_width/3))

			if ((i+1)/3) == np.floor((i+1)/3) and i>0:
				binned_signal.append(binsize[int(i-2)] + binsize[int(i-1)] + binsize[int(i)])
				sumbins = sumbins + binned_signal[int(i/3)]
				startingfreqs.append(scanfreq-bin_width)

			i+=1
			if len(binned_signal)>500:
				break

	elif shape_model == "pCDM_only":
		modf = float(((vsol+vel_mod)**2)/(vsol**2))
		while sumbins<cutoffarea:
			scanfreq = float(startfreq+(i)*(bin_width/3))*10**(-6)
			modulatingfreq = ((scanfreq - RMF)*modf+RMF)
			binsize.append(pCDM_only(m, modulatingfreq, vsol)*(bin_width/3))

			if ((i+1)/3) == np.floor((i+1)/3) and i>0:
				binned_signal.append(binsize[int(i-2)] + binsize[int(i-1)] + binsize[int(i)])
				sumbins = sumbins + binned_signal[int(i/3)]
				startingfreqs.append(scanfreq-bin_width)
			i += 1
			if len(binned_signal)>500:
				break

	elif shape_model == "pCDM_w_baryons":
		modf = ((vsol+vel_mod)**2)/(vsol**2)
		while sumbins<cutoffarea:
			scanfreq = float(startfreq+(i)*(bin_width/3))*10**(-6)
			modulatingfreq = ((scanfreq - RMF)*modf+RMF)
			binsize.append(pCDM_w_baryons(m, modulatingfreq, vsol)*(bin_width/3))
			if ((i+1)/3) == np.floor((i+1)/3) and i>0:
				binned_signal.append(binsize[int(i-2)] + binsize[int(i-1)] + binsize[int(i)])
				sumbins = sumbins + binned_signal[int(i/3)]
				startingfreqs.append(scanfreq-bin_width)
			i += 1
			if len(binned_signal)>500:
				break

	elif shape_model == "pCDM_maxwell_like_form":
		modf = ((vsol+vel_mod)**2)/(vsol**2)
		while sumbins<cutoffarea:
			scanfreq = float(startfreq+(i)*(bin_width/3))*10**(-6)
			modulatingfreq = ((scanfreq - RMF)*modf+RMF)
			binsize.append(pCDM_maxwell_like_form(m, modulatingfreq, vel_mod, alpha, beta, T)*(bin_width/3))

			if ((i+1)/3) == np.floor((i+1)/3) and i>0:
				binned_signal.append(binsize[int(i-2)] + binsize[int(i-1)] + binsize[int(i)])
				sumbins = sumbins + binned_signal[int(i/3)]
				startingfreqs.append(scanfreq-bin_width)
			i += 1
			if len(binned_signal)>500:
				break

	elif shape_model == "axionDM_only":
		while sumbins<cutoffarea:
			modf = ((vsol+vel_mod)**2)/(vsol**2)
			scanfreq = float(startfreq+(i)*(bin_width/3))*10**(-6)
			modulatingfreq + ((scanfreq - RMF)*modf+RMF)
			binsize.append(axionDM_only(m, modulatingfreq)*(bin_width/3))

			if ((i+1)/3) == np.floor((i+1)/3) and i>0:
				binned_signal.append(binsize[int(i-2)] + binsize[int(i-1)] + binsize[int(i)])
				sumbins = sumbins + binned_signal[int(i/3)]
				startingfreqs.append(scanfreq-bin_width)
			i += 1
			if len(binned_signal)>500:
				break

	#Error message if shape model has not yet been defined or doesnt even exist
	if shape_model != "axionDM_only" and shape_model != "pCDM_maxwell_like_form" and shape_model != "pCDM_w_baryons" and shape_model != "pCDM_only" and shape_model != "SHM" and shape_model != "axionDM_w_baryons":
		return "Error: Axion shape model not recognized"

	#if wantTseries == 'y', then fourier transform the signal into a time series. if not, leave it alone
	if wantTseries == "y":
		finalsignal = np.fft.fft(binned_signal)
	elif wantTseries == "n":
		finalsignal = binned_signal

	#collect important properties and values
	info = {'rest mass frequency':RMF, 'shape': shape_model, 'signal': finalsignal, 'freqs':startingfreqs}

	return info


"""
plt.style.use('seaborn-pastel')
alpha = 0.5
lgndfntsz = 11
lw = 3.0
fntsz = 18

plotall = input('Do you want to plot a signal?: (y/n)')
 #Non-functioning...

if plotall == 'y':
	shape_model = input('What do you want to plot?: ("SHM","axionDM_w_baryons", "pCDM_only", "pCDM_w_baryons", "pCDM_maxwell_like_form", "axionDM_only"):   ')
	mod_type = input('What modulation do you want to observe?: ("solar", "earth orbit", "earth rotation")   ')
	timestamp = repr(dt.strftime(dt.now(), "%Y-%m-%d %H:%M:%S-%f"))
	m = 3*10**(-6)
	res = 100
	wantTseries = 'n'
	startfreq = None
	alpha = 0.3
	beta = 1.3
	T = 4.7*10**(-7)

	info = modulatedsignal(mod_type, timestamp, shape_model, m, res, wantTseries, startfreq, alpha, beta, T)
	c = Cnsts.c.value*3*10**(-3)
	RME = 10**(-6)
	m = RME/(c**2)
	vsol = 232
	h = Cnsts.h.value*6.242*10**18
	fo = RME/h
	domain = info['freqs']
	signal = info['signal']

	plt.figure()
	plt.xlabel('MHz')
	plt.ylabel('Probability density')
	plt.xlim(np.min(domain), np.max(domain))
	plt.ylim(0, np.max(signal)+(np.max(signal)/5))
	plt.plot(domain, signal, linestyle=":", alpha = alpha, color='green')
	plt.legend(loc = 0, shadow=True)
	plt.show()

	print(shape_model, mod_type, timestamp)
"""
"""
m = 3*10**(-6)
h = Cnsts.h.value*6.242*10**18
rmf = m/h*10**(-6)
f = (m/h + 125)*10**(-6)
print(modulatedsignal('solar', '2018-6-21 14:26:21-22', 'SHM', m, 100,None,None,0.3,1.3,4.7*10**(-7)))
"""
