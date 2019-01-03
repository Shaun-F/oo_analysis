"""
modulation.py: methods for bringing axion spectra into the lab reference frame

Created by: Erik Lentz
Creation Date: 6/1/18
"""
import numpy as np
import scipy as sp
from astropy import constants as Cnsts
from datetime import datetime as dt
from dateutil.parser import parse
import pytz
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.pylab as plt
import pandas as pd
import argparse

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
		modulation_type = self.pec_vel
		timestamps = self.timestamp
		shape_model = self.signal
		axion_mass = self.axion_mass
		scans = self.scans
		
		signals = {}
		for key in scans.keys():
			signals[key] = self.modulatedsignal(modulation_type, timestamps[key], shape_model, axion_mass)
		return signals
	# methods for incorporating several levels of peculiar velocities and other
	#modulating effects
	def Equ_to_Gal(self,numberofhours):
		"""
		Description: Function returns the velocity vector of the lab in galactic 	cylindrical coordinates using a matrix equation found in Spherical Astronomy by Robin Green
		Parameters: Number of hours since B1950
		Output:Velocity vector in Galactic cylindrical Coordinates
		"""

		phi = 123 # angle between vr and vt at the start of B1950 epoch
		v=0.4638*np.cos(47.66) #0.4638 is rotation speed of earth, 47.66 is latitude of experiment

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

	def earth_orbit_velocity(self, date):
		"""
		Description: Function returns the orbital velocity in the cylindrical galactic center basis vr, vt, vz
		Parameters:Timestamp
		Output: Dictionary of velocity components 'vr', 'vt', 'vz'
		"""
		time = timestamptranslator(date)
		t = time*86164.09056/86400
		T = 3.1558*10**7
		phi = 86.14
		v = 29.785
		v_plane = {"vx":v*np.cos(2*np.pi*t/(T+phi)), "vy":v*np.sin(2*np.pi*t/(T+phi)), "vz":0}
		v_gal = {"vr":Equ_to_Gal(time)["vr"]*v_plane["vx"], "vt":Equ_to_Gal(time)["vt"]*v_plane["vy"]}

		return v_gal

	def cavity_rotation_velocity(self, date):
		"""
		Description: Function returns the cavitys orbital velocity in cylindrical galactic center basis vr, vt, vz
		Parameters: timestamp
		Output: Dictionary to velocity components 'vr', 'vt', 'vz'
		"""
		time = timestamptranslator(date)
		t = time*3600
		T = 86164.09056
		phi = -122.3 #phi is the angle between vr and vt  on June 17th 1998 -- because the earth is closest to the center of galaxy on this day
		v = 0.4638*np.cos(47.66)
		v_spin = {'vx': v*np.cos(2*np.pi*t/(T+phi)), 'vy':v*np.cos(2*np.pi*t/(T+phi)), 'vz':0}
		v_gal = {'vy': Equ_to_Gal(time)['vr']*v_spin['vx'], 'vt': Equ_to_Gal(time)['vt']*v_spin['vy']}

		return v_gal

	def modulate(self,modulation_type, date):
		hourtime = timestamptranslator(timestamp)

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

	def modulatedsignal(self,modulation_type, timestamp, shape_model, axion_mass, resolution=None, wantTseries=None, startfreq=None, alpha=None, beta=None, T=None):
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
		vel_mod = float(self.modulate(modtype, timestamp)) #variation of experimental speed around mean solar speed
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
