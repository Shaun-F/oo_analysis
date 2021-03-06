"""
modulation.py: methods for bringing axion spectra into the lab reference frame

Created by: Erik Lentz
Creation Date: 6/1/18
"""
import numpy as np; import numpy
from astropy import constants as Cnsts
from astropy import units as U
from datetime import datetime as dt
from dateutil.parser import parse
import pytz
import sys
import oo_analysis.toolbox.coord_transform as tf_lib
from oo_analysis.signals.signal_lib import signal
from oo_analysis.toolbox.signal_width import calc_signal_width

class modulation():
	def __init__(self,**params):
		for key,param in params.items(): #python 3.x
			setattr(self,key,param)
		
		return None
	def executor(self):
		"""
		Function runs modulation routine 
		"""
		#Define default values
		default_keys = ["pec_vel", "timestamp", "signal", "axion_mass", "dig_dataset"]
		default = {key: getattr(self, key) for key in default_keys}
		
		#Pull scan properties
		modulation_type = default["pec_vel"]
		timestamps = default["timestamp"]
		shape_model = default["signal"]
		axion_mass = default["axion_mass"]
		scans = default["dig_dataset"]
		
		self.mod_vels = self.modulate(modulation_type, timestamps, shape_model) #dictionary of total experiment velocity around galactic center keyed by timestamp
		

		#define secondary values and attach to class
		secondary={}
		for key in self.__dict__.keys():
			if key not in default_keys and key!='mod_vels':
				secondary[key] = getattr(self, key)
				
		signals = {}
		if isinstance(scans, list) or isinstance(scans, numpy.ndarray) or isinstance(scans, dict): #If more than one scan is being analyzed
			counter = 1
			N_iter = len(self.keys)
			for key in self.keys:
				if not scans[key].attrs['cut']: #No need to generate signals for scans that wont be analyzed
					timestamp = timestamps[key]
					if not isinstance(axion_mass, float):
						a_mass = axion_mass[key]
					else:
						a_mass = axion_mass
					signals[key] = self.modulatedsignal(timestamp, shape_model, a_mass, self.mod_vels[timestamp], **secondary) #Generate modulated signals
					if ('SIA',True) not in list(self.__dict__.items()):
						if int(counter/N_iter*100)-counter/N_iter*100<10**(-8):
							print("generating signals ({0} % complete)                 \r".format(int((counter/N_iter)*100)), end='')
					counter+=1
				else:
					pass
		else: #If only a single scan is being analyzed
			if not scans.attrs['cut']:
				scan_id = scans.name[-6:] #Name usually includes the directory too, so just splice out the 6-digit id
				timestamp = timestamps[scan_id]
				a_mass = axion_mass
				signals[scan_id] = self.modulatedsignal(timestamp, shape_model, a_mass, self.mod_vels[timestamp], **secondary)
			else:
				pass
		return signals
		
	# methods for incorporating several levels of peculiar velocities and other
	#modulating effects
	
	def modulate(self,modulation_type, timestamp, shape_model):

		#vt_orbit = vel_orbit["vt"]
		#vt_rotation = vel_rotation["vt"]

		if modulation_type == None or modulation_type =="None":
			modulation_type = "solar"
		
		if shape_model!="axionDM_w_baryons": #N-Body signal shape requires the peculiar velocity of earth around sun, not earths total velocity
			if modulation_type == "solar":
				return tf_lib.transformer().solar_vel_GalacticFrame(timestamp)
			elif modulation_type=="earth orbit":
				return tf_lib.transformer().earth_vel_GalacticFrame(timestamp)
			elif modulation_type=="earth rotation":
				return tf_lib.transformer().experiment_vel_GalacticFrame(timestamp)
			else:
				return "Error: modulation type not recognized"
		else:
			solar = tf_lib.transformer().solar_vel_GalacticFrame(timestamp)
			if modulation_type=="solar":
				return {time: [0,0,0] for time in list(timestamp.values())}
			elif modulation_type=='earth orbit':
				earth_vel = tf_lib.transformer().earth_vel_GalacticFrame(timestamp)
				return {key: np.array(earth_vel[key]) - np.array(solar.get(key,0)) for key in list(solar.keys())}
			elif modulation_type=='earth rotation':
				experiment_vel = tf_lib.transformer().experiment_vel_GalacticFrame(timestamp)
				return {key: np.array(experiment_vel[key]) - np.array(solar.get(key,0)) for key in list(solar.keys())}

	def modulatedsignal(self, timestamp, shape_model, axion_mass, vel_mod, **kwargs):
		"""
		Description:Function bins a given velocity-modulated axion shape by integration.

		Parameters: 
			Modulation types include 'solar', 'earth orbit', 'earth rotation'. 
			Timestamp takes either the form "dd-mm-yyyy  hh:mm:ss AM", which is used by the Digitizer, or the form "yyyy-mm-dd hh:mm:ss-ss", which is the standard ISO 8601 format. 
			Current shape models include "SHM", "axionDM_w_baryons", "pCDM_only", "pCDM_w_baryons", "pCDM_maxwell_like_form", "axionDM_only". 
			Axion mass is in eV.

		Optional Parameters: 
			Resolution give the size of the bins in Hz. 
			wantTseries specifies if the output should be in frequency space or as a timeseries (y/n). 
			startfreq allows you to specify the position of the signal shape within the array (MHz). 
			alpha, beta, T are fit parameters for the N-Body model.

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
		if 'mod_vels' in kwargs.keys():
			vel_mod = kwargs['mod_vels']

		#modulation velocity
		#specify parameters to be used in various signal scripts
		vsol = tf_lib.transformer().solar_vel_GalacticFrame() # (km/s) mean solar speed around galaxy
		h = Cnsts.h.to(U.eV*U.s).value # plancks constant in eV*s
		m = axion_mass # axion mass in eV
		bin_width = float(resolution) #size of bins in hZ
		RMF = float(m/h) # Rest mass frequency in hZ
		cutoffarea = 0.92
		

		signal_cl = signal(**self.__dict__)

		#set default value for starting frequency
		if 'startfreq' not in kwargs.keys():
			startfreq = RMF

		binsize = []
		binned_signal = []
		sumbins = 0

		startingfreqs = [startfreq]
		info = []
		i=0

		#modulate specified axion shape based on input parameters
		if shape_model == "axionDM_w_baryons":
			while sumbins<cutoffarea:
				scanfreq = float(startfreq + (i)*(bin_width/3))*10**(-6)
				binsize.append((signal_cl.axionDM_w_baryons(scanfreq, mod_vel = vel_mod, v_sol = vsol, mass = m)*(bin_width/3))[0])
				if ((i+1)/3) == np.floor((i+1)/3) and i>0:
					binned_signal.append(binsize[int(i-2)] + binsize[int(i-1)] + binsize[int(i)])
					sumbins = sumbins + binned_signal[int(i/3)]
					startingfreqs.append(scanfreq*10**(6)-bin_width)
				i+=1
				if len(binned_signal)>500:
					break

		elif shape_model == "SHM":
			while sumbins<cutoffarea:
				scanfreq = float(startfreq+(i)*(bin_width/3))*10**(-6)
				binsize.append((signal_cl.SHM(scanfreq, m, vel_mod)*(bin_width/3))[0])

				if ((i+1)/3) == np.floor((i+1)/3) and i>0:
					binned_signal.append(binsize[int(i-2)] + binsize[int(i-1)] + binsize[int(i)])
					sumbins = sumbins + binned_signal[int(i/3)]
					startingfreqs.append(scanfreq*10**(6)-bin_width)

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
					startingfreqs.append(scanfreq*10**(6)-bin_width)
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
					startingfreqs.append(scanfreq*10**(6)-bin_width)
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
					startingfreqs.append(scanfreq*10**(6)-bin_width)
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
					startingfreqs.append(scanfreq*10**(6)-bin_width)
				i += 1
				if len(binned_signal)>500:
					break
		else:
			return "Error: Axion shape model not recognized" #Error message if shape model has not yet been defined or doesnt even exist

		
		#if wantTseries == 'y', then fourier transform the signal into a time series. if not, leave it alone
		if wantTseries == "y":
			finalsignal = np.fft.fft(binned_signal)
		elif wantTseries == "n":
			finalsignal = binned_signal

		#collect important properties and values
		if len(finalsignal)>100: #upper bound on length of signal array to mitigate floating point error effects
			finalsignal = finalsignal[:100]
			
		try:
			signal_width = calc_signal_width(np.array(finalsignal))*bin_width #This is in hertz!
		except IndexError:
			print(finalsignal)
		info = {'rest mass frequency':RMF, 
				'shape': shape_model, 
				'signal': np.asarray(finalsignal), 
				'freqs':np.asarray(startingfreqs),
				'signal_width': signal_width}
				
		return info

	def timestamptranslator(self, date):
		"""
		Description:Function takes in a time stamp either in the form 'dd-mm-yyyy  hh:mm:ss AM' or 'yyyy-mm-dd hh:mm:ss-ss' and returns the total hours elapsed since June 17th, 1998 midnight GMT (Input timezone is PST).
		Parameters: timestamp
		Output: Number of hours since Midnight june 17th, 1998 GMT from input timestamp

		"""

		
		d1 = parse(date)
		d1 = dt(d1.year,d1.month,d1.day,d1.hour,d1.minute,d1.second,0, tzinfo=pytz.timezone('America/Vancouver'))
		tz_pst = pytz.timezone('America/Vancouver') #Pacific standard time
		tz_grnwch = pytz.timezone('Etc/Greenwich') #greenwich mean time

		d2 = dt(1949,12,31,23,58,56,3795, tzinfo = tz_grnwch) #B1950 epoch in UTC, the standard UNIX reference date
		return ((d1-d2).total_seconds()/(60*60)) #Number of hours since B1950
