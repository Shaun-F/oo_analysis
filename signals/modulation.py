"""
modulation.py: methods for bringing axion spectra into the lab reference frame

Created by: Erik Lentz
Creation Date: 6/1/18
"""
import numpy as np
from astropy import constants as Cnsts
from datetime import datetime as dt
from dateutil.parser import parse
import pytz
import sys
sys.path.append("..")
import toolbox.coord_transform as tf_lib
from signals.signal_lib import signal

class modulation():
	def __init__(self,**params):
		for key,param in params.items(): #python 3.x
			setattr(self,key,param)
		# sets up time, peculiar velocities, etc
		return None
	def executor(self):
		"""
		Function runs modulation routine 
		"""
		#Define default values
		default_keys = ["pec_vel", "timestamp", "signal", "axion_mass", "dig_dataset"]
		default = {key: getattr(self, key) for key in default_keys}
		
		#define secondary values
		secondary={}
		for key in self.__dict__.keys():
			if key not in default_keys:
				secondary[key] = getattr(self, key)
			
		modulation_type = default["pec_vel"]
		timestamps = default["timestamp"]
		shape_model = default["signal"]
		axion_mass = default["axion_mass"]
		scans = default["dig_dataset"]
		
		signals = {}
		for key in self.keys:
			signals[key] = self.modulatedsignal(modulation_type, timestamps[key], shape_model, axion_mass, **secondary)
		return signals
		
	# methods for incorporating several levels of peculiar velocities and other
	#modulating effects
	
	def modulate(self,modulation_type, timestamp):

		vel_orbit = tf_lib.transformer().earth_vel_GalacticFrame(timestamp)
		#vt_orbit = vel_orbit["vt"]

		vel_rotation = tf_lib.transformer().experiment_vel_GalacticFrame(timestamp)
		#vt_rotation = vel_rotation["vt"]

		vel_solar = tf_lib.transformer().solar_vel_GalacticFrame()
		if modulation_type == None or modulation_type =="None":
			modulation_type = "solar"

		if modulation_type == "solar":
			return vel_solar
		elif modulation_type=="earth orbit":
			return vel_orbit
		elif modulation_type=="earth rotation":
			return vel_rotation 
		else:
			return "Error: modulation type not recognized"

	def modulatedsignal(self, modulation_type, timestamp, shape_model, axion_mass, **kwargs):
		"""
		Description:Function bins a given velocity-modulated axion shape by integration.

		Parameters: Modulation types include 'solar', 'earth orbit', 'earth rotation'. Timestamp takes either the form "dd-mm-yyyy  hh:mm:ss AM", which is used by the Digitizer, or the form "yyyy-mm-dd hh:mm:ss-ss", which is the standard ISO 8601 format. Current shape models include "SHM", "axionDM_w_baryons", "pCDM_only", "pCDM_w_baryons", "pCDM_maxwell_like_form", "axionDM_only". Axion mass is in eV.

		Optional Parameters: Resolution give the size of the bins in Hz. wantTseries specifies if the output should be in frequency space or as a timeseries (y/n). startfreq allows you to specify the position of the signal shape within the array (MHz). alpha, beta, T are fit parameters for the N-Body model.

		Output: velocity-modulated signal shape over frequency or time.
		"""

		# Set default values for signals not being used
		if 'resolution' not in kwargs.keys():
			resolution = 100
		else: 
			resolution = kwargs["resolution"]
		if 'wantTseries' not in kwargs.keys():
			wantTseries='n'
		else:
			wantTseries = kwargs["wantTseries"]
		if 'alpha' not in kwargs.keys():
			alpha = 0
		else:
			alpha = kwargs["alpha"]
		if 'beta' not in kwargs.keys():
			beta = 0
		else:
			beta = kwargs["beta"]
		if 'T' not in kwargs.keys():
			T = 0
		else:
			T = kwargs["T"]

		#specify parameters to be used in various signal scripts
		modtype = modulation_type #create pointer with shorter character length
		vsol = tf_lib.transformer().solar_vel_GalacticFrame() # (km/s) mean solar speed around galaxy
		vel_mod = float(self.modulate(modtype, timestamp)) #variation of experimental apparatus speed around sun w.r.t. galactic center
		h = Cnsts.h.value*6.242*10**18 # plancks constant in eV*s
		m = axion_mass # axion mass in eV
		bin_width = float(resolution) #size of bins in hZ
		RMF = float(m/h) # Rest mass frequency in hZ
		cutoffarea = 0.9

		signal_cl = signal(self.__dict__)

		#set default value for starting frequency
		if 'startfreq' not in kwargs.keys():
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
				binsize.append(signal_cl.axionDM_w_baryons(vel_mod, vsol, m, scanfreq)*(bin_width/3))
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
				binsize.append(signal_cl.SHM(scanfreq, m, vel_mod)*(bin_width/3))

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
				binsize.append(signal_cl.pCDM_only(m, modulatingfreq, vsol)*(bin_width/3))

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
				binsize.append(signal_cl.pCDM_w_baryons(m, modulatingfreq, vsol)*(bin_width/3))
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
				binsize.append(signal_cl.pCDM_maxwell_like_form(m, modulatingfreq, vel_mod, alpha, beta, T)*(bin_width/3))

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
				binsize.append(signal_cl.axionDM_only(m, modulatingfreq)*(bin_width/3))

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

	def timestamptranslator(self, date):
		"""
		Description:Function takes in a time stamp either in the form 'dd-mm-yyyy  hh:mm:ss AM' or 'yyyy-mm-dd hh:mm:ss-ss' and returns the total hours elapsed since June 17th, 1998 midnight GMT (Input timezone is PST).
		Parameters: timestamp
		Output: Number of hours since Midnight june 17th, 1998 GMT from input timestamp

		"""
		
#		if self.timestamp[2]== '-':
#			"""This IF statement tests whether the date string is in the form 'dd-mm-yyyy  hh:mm:ss AM'. This is the case for a digitizer timestamp"""
#			d1 = dt.strptime(timestamp, "%d-%m-%Y  %I:%M:%S %p")
#			d1 = dt(d1.year,d1.month,d1.day,d1.hour,d1.minute,d1.second,0, tzinfo=pytz.timezone('America/Vancouver'))
#		if self.timestamp[4]== '-':
#			"""This IF statement tests whether the date string is in the form 'yyyy-mm-dd hh:mm:ss-ss'. This is the standard ISO 8601 format and the format Python uses for the datetime module"""
#			d1 = dt.strptime(timestamp, "%Y-%m-%d %H:%M:%S-%f")
#			d1 = dt(d1.year,d1.month,d1.day,d1.hour,d1.minute,d1.second,0, tzinfo=pytz.timezone('America/Vancouver'))
		
		d1 = parse(date)
		d1 = dt(d1.year,d1.month,d1.day,d1.hour,d1.minute,d1.second,0, tzinfo=pytz.timezone('America/Vancouver'))
		tz_pst = pytz.timezone('America/Vancouver') #Pacific standard time
		tz_grnwch = pytz.timezone('Etc/Greenwich') #greenwich mean time

		d2 = dt(1949,12,31,23,58,56,3795, tzinfo = tz_grnwch) #B1950 epoch in UTC, the standard UNIX reference date
		return ((d1-d2).total_seconds()/(60*60)) #Number of hours since B1950
