"""
filter_lib.py: holds axion filter types

Created by: Erik Lentz
Creation Date: 6/1/18
"""
import numpy as np
import scipy as sp
from astropy import constants as Cnsts
from astropy import units as U
from datetime import datetime as dt
import pytz
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.pylab as plt
import pandas as pd
import argparse


class signal(object):
	# want to be able to use as part of the mainstream analysis and as a library
	def __init__(self,params):
		for key,param in params.items(): #python 3.x
			setattr(self,key,param)
		#callable = getattr(self,self.signal) # or whatever the rel. attr. is
		return None


	# series of defs giving base axion spectra
	def top_hat(self):
		# number of bins wide sets the distribution
		n_bins = self.nbins
		return np.ones(n_bins)/n_bins

	def SHM(self,freq, mass, modulation_vel):
		"""
		Modulation_vel is the total velocity of the experiment around the galactic center
		"""
		sigma_v = 160
		v_lab = modulation_vel

		return self.ITS(freq, mass, sigma_v, v_lab)

	def pCDM_only(self,mass,freq, v_lab):
		"""
		Description: Function returns power spectra for pressureless cold dark matter without baryons
		Parameters:mass of axion (eV), Sampling frequency (Hz)
		Output: Power spectra at sampling frequency
		"""
		alpha = 0.25
		beta = 1.3
		T = 2.0*10**(-7)

		return pCDM_maxwell_like_form(mass,freq, v_lab, alpha,beta,T)

	def pCDM_w_baryons(self,mass,freq, v_lab):
		"""
		Description: Function returns power spectra for pressureless cold dark matter, taking into account baryons
		Parameters: Mass of axion (eV), Sampling frequency (Hz)
		Output: Power spectra at sampling frequency
		"""
		alpha = 0.36
		beta = 1.39
		T = 4.7*10**(-7)

		return pCDM_maxwell_like_form(mass,freq,v_lab,alpha,beta,T)

	def pCDM_maxwell_like_form(self, mass, freq, v_lab, vsol, alpha, beta, T):
		"""
		Description:Function gives the frequency power spectra for halo axions as modeled by the N-Body simulation using pCDM, parameterized by three values alpha, beta, T
		Parameters:mass of axion (eV), sampling frequency (Hz), velocity of laboratory (vector in km/s), solar velocity (vector in km/s), parameter values of N-Body simulation
		Output: Power spectra of axion
		"""
		f = freq*10**6 #sampling frequency
		c = Cnsts.c.value*10**(-3) #Speed of light
		RME = mass #Rest mass energy of axion
		h = Cnsts.h.value*6.242*10**18 # plancks constant in ev*s
		fo = RME/h # Rest mass frequency of axion
		
		v_pec = v_lab-vsol #peculiar velocity of experiment around sun
		gamma = sp.special.gamma((1.0+alpha)/beta)
		Cnum = beta
		Boost = 1-sum(v_pec*vsol)/((vsol**2).sum())
		Cden = ((Boost/(T*fo))**beta)**(-(1.0+alpha)/beta)*(Boost/(T*fo))**(alpha)*gamma
		Cnst = Cnum/Cden

		if fo<f:
			A = Cnst*(h*(f-fo)/(mass*T))**(alpha)*np.exp(-(h*(f-fo)/(mass*T))**beta)
		else:
			A = 0

		return A

	def ITS(self, freq, mass, sigma_v, v_lab): #bin_width,sigmav,vt,vr, freq,mass,sigma_v,v_lab
		"""Description: Function takes in parameters and outputs probability to find axion at input frequency
		Parameters: frequency (MHz) to find corresponding probability. Mass of axion in ev. velocity dispersion of axions in Kilometers per second. velocity of laboratory in galactic frame in Kilometers per second.
		Output: Probability to find axion at input frequency"""
		if isinstance(freq, float) or isinstance(freq, int):
			freq = [freq] #convert to list
		freq_arr = np.asarray(freq) * 10**(6) #in Hz
		
		c = Cnsts.c.to(U.km/U.s).value #Speed of light in kilometers per second
		h = Cnsts.h.to(U.eV*U.s).value #Plancks constant in units eV seconds
		RME = mass #Rest Mass energy of Axion in eV
		rmfreq = (RME/h) # Rest mass frequency of axion
		v = np.sqrt(2) * c * np.sqrt((freq_arr*h)/(RME) - 1) #velocity of axion
		v_lab_mag = np.linalg.norm(v_lab) #speed of the experiment in a right handed coordinate system center on the galactic center
		beta = 1/(2*(sigma_v**2)) #this is the beta from turner 1990
		mask = np.where(freq_arr>rmfreq)
		mask_complement = np.setdiff1d(range(len(freq_arr)), mask) #If input frequency is less than the rest mass frequency, probability is zero (nonphysical)
		v_store = v.copy()

		v_store[mask] = (2*h*(c**2)*np.sqrt(beta/np.pi))/(RME*v_lab_mag) * np.exp(-beta*(v_lab_mag**2+v[mask]**2)) * np.sinh(2*beta*v_lab_mag*v[mask])
		v_store[mask_complement] = 0 #If input frequency is less than the rest mass frequency, probability is zero (nonphysical)
		return v_store

	def axionDM_w_baryons(self, mod_vel, v_sol, mass, freq):

		"""
		Description: Function gives probability density for axion signal of Lentz et al. 2017
		Parameters: mod_vel (km/s) is the variation of the experiments velocity from the mean solar velocity w.r.t a right-handed galactic coordinate system centered on the galactic center. v_sol (km/s) is the mean solar velocity. mass (eV) is the mass of the axion. freq (MHz) is a desired sample point.
		Output: Probability to find an axion with given input frequency.
		"""
		
		#NO need for the following. velocities are now vectored quantities
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
		"""

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

		Boost = 1-sum((mod_vel)*(v_sol))/((v_sol**2).sum()) #Modulating parameter to account for diurnal variations



		N_num = beta
		N_den = ((Boost/(T*f_o))**beta)**(-(1.0+alpha)/beta)*(Boost/(T*f_o))**(alpha)*gamma
		nrmcnst = N_num/N_den

		#The distribution function
		A = nrmcnst*(((f-f_o)*Boost)/(f_o*T))**(alpha)*np.exp(-(((f-f_o)*Boost)/(f_o*T))**beta)

		return A
