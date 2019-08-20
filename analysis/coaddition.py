
import datetime as dt
import numpy as np
import sys
import os

from oo_analysis.toolbox.add_to_dataset import * #addtodataset, subtractfromdataset, assign_newdata
import numba; from numba import njit, jit
from oo_analysis.toolbox.grand_spectra_init import initialization


class scan_cl(object):
	"""
	This class performs the necessary addition or subtraction of scans from a grand spectrum.
	Attributes: scan_number, ch, res, modulation_type, axion_shape
	Methods: axion_scan, digitizer_scan, analyze_scan, add_scan
	"""
	
	
	def __init__(self, scid, scan, chunk, op, **params):
		"""
		Initialization method
		"""
		
		self.op = op #Whether to add or subtract scan from grand spectra
		self.scan=scan #Scan dictionary
		self.scan_scan_number = scid #Scan number
		
		#setup class definitions for chunk 
		#pull only array of points from grand_spectra that line up with scan
		self.chunk = chunk #Grand spectra HDF5 dataset
	
	def axion_fit_consolidation(self, inx_start, inx_end):
		"""
		Coadd scan axion fit into grand spectra
		"""
		#return divider(self, self.chunk['optimal_weight_sum'][inx_start:inx_end], self.chunk['model_excess_sqrd'][inx_start:inx_end])

		return inverse_quadrature(self, self.chunk['axion_fit'][inx_start:inx_end], self.scan['axion_fit'], self.op, -2)
		
	def axion_fit_significance_consolidation(self, inx_start, inx_end):
		"""
		Coadd scan axion fit significance  into grand spectra. Since axion fit and axion fit uncertainty are already coadded, can simply divide the grand spectra axion fit and axion fit 
		uncertainty to get the coadded axion fit significance
		"""
		
		return divider(self, self.chunk['axion_fit'][inx_start:inx_end], self.chunk['axion_fit_uncertainty'][inx_start:inx_end])
	
	def coupling_sensitivity_consolidation(self, inx_start, inx_end):
		"""
		Coadd the scan coupling sensitivity into the grand spectra. 
		"""
		
		return inverse_root_quadrature(self, self.chunk['sensitivity_coupling'][inx_start:inx_end], self.scan['sensitivity_coupling'], self.op, -4)
	
	def maximum_likelihood_uncertainty_consolidation(self, inx_start, inx_end):
		"""
		Coadded the scan maximum likelihood uncertainty into the grand spectra
		"""
		
		return inverse_quadrature(self, self.chunk['axion_fit_uncertainty'][inx_start:inx_end], self.scan['axion_fit_uncertainty'], self.op, -2)
	
	def model_excess_sqrd_consolidation(self, inx_start, inx_end):
		"""
		Coadd the scans model excess sqrd (quadratic coefficient of chi squared) into the grand spectra 
		"""
		
		return add_sub(self, self.chunk['model_excess_sqrd'][inx_start:inx_end], self.scan['model_excess_sqrd'], self.op)
	
	def noise_power_consolidation(self, inx_start, inx_end):
		"""
		Coadd the scans noise power into the grand spectra. Since SNR and the optimal weight sum were already coadded, can just us the grand spectra SNR and optimal weight sum
		"""
		
		return divider(self, self.chunk['SNR'][inx_start:inx_end], self.chunk['optimal_weight_sum'][inx_start:inx_end])
	
	def nscan_consolidation(self, inx_start, inx_end):
		"""
		Coadd the scans nscan into the grand spectra
		"""
		
		return add_sub(self, self.chunk['nscans'][inx_start:inx_end], self.scan['nscans'], self.op)
		
	def optimal_weight_sum_consolidation(self, inx_start, inx_end):
		"""
		Coadd the scans optimal weight sum (linear coefficient of chi squared) into the grand spectra
		"""
		
		return add_sub(self, self.chunk['optimal_weight_sum'][inx_start:inx_end], self.scan['optimal_weight_sum'], self.op)
	
	def power_deviation_consolidation(self, inx_start, inx_end):
		"""
		Coadd the scans power deviations into the grand spectra
		"""
		
		return add_sub(self, self.chunk['power_deviation'][inx_start:inx_end], self.scan['power_deviation'], self.op)
	
	def power_sensitivity_consolidation(self, inx_start, inx_end):
		"""
		Coadd the scans power sensitivity into the grand spectra
		"""
		
		return inverse_quadrature(self, self.chunk['sensitivity_power'][inx_start:inx_end], self.scan['sensitivity_power'], self.op, -2)
		
	def sigma_A_consolidation(self, inx_start, inx_end):
		"""
		Coadd the scans axion fit uncertainty into the grand spectra
		"""
		
		return inverse_quadrature(self, self.chunk['axion_fit_uncertainty'][inx_start:inx_end], self.scan['axion_fit_uncertainty'], self.op, -2)
	
	def SNR_consolidation(self, inx_start, inx_end):
		"""
		Coadd the scans SNR into the grand spectra
		"""
		
		return quadrature(self, self.chunk['SNR'][inx_start:inx_end], self.scan['SNR'], self.op, 2)

	"""
	def weighted_delta_consolidation(self):
		
		weighted_deltas = self.chunk_weighted_deltas
		swds = self.scan_weighted_deltas
		for i, val in enumerate(weighted_deltas):

			if self.op=="+":
				weighted_deltas[i] = val + swds[i]
			else:
				weighted_deltas[i] = val - swds[i]
		
		return add_sub(self, self.chunk_weighted_deltas, self.scan['weighted_deltas'], self.op)
	"""
	def close_out(self, chunk):
		"""
		DEPRICATED. I coadd directly to the HDF5 datasets instead.
		Class method saves the coadded data to the grand spectra. 
		"""
		chunk['nscans'][...] = self.chunk['nscans']
		chunk['sigma_w'][...] = self.chunk['sigma_w']
		chunk['optimal_weight_sum'][...] = self.chunk['optimal_weight_sum']
		chunk['SNR'][...] = self.chunk['SNR']
		chunk['noise_power'][...] = self.chunk['noise_power']
		chunk['power_deviation'][...] = self.chunk['power_deviation']
		chunk['model_excess_sqrd'][...] = self.chunk['model_excess_sqrd']
		chunk['axion_fit'][...] = self.chunk['axion_fit']
		chunk['axion_fit_uncertainty'][...] = self.chunk['axion_fit_uncertainty']
		chunk['sensitivity_power'][...] = self.chunk['sensitivity_power']
		chunk['sensitivity_coupling'][...] = self.chunk['sensitivity_coupling']
		chunk['last_change'][...] = self.chunk['last_change']
		print(self.chunk['last_change'])

		

def add_subtract_scan(add_subtract, scan, object, scan_id, grand_spectra_group, **params):
	"""
	Parameters
		add_subtract: ('add', 'subtract') Determines what to do with scan
		scan: (dictionary) items to add to grand spectra, including deltas
		object: (class) coaddition class thats adds scan into grand spectra
		scan_id: (int) id of scan
		grand spectra group: (HDF5 group) the grand spectra group to be coadded into
	"""
	#initialize scan object	
	if add_subtract == "add":
		object.op = "+"
	elif add_subtract == "subtract":
		object.op = "-"
	else:
		return "Error: combining operation not recognized. Available operations are 'add' and 'subtract'"
	
	
	corrupt_scan=False
	if type(scan) is str:
		corrupt_scan = True
		add_subtract = 'ommit'
		object.op = 'nil'
	
	#update class attributes
	object.scan = scan
	object.scan_scan_number = scan_id
	
	#Get indices of grand spectra where scan will be added into
	appended_length = len(scan['axion_frequencies'])
	mid_freq = scan['middle_frequency']
	ainx_start = int(numpy.round((mid_freq - grand_spectra_group['axion_frequencies'][0])/(95.4)) - appended_length/2)
	ainx_end = int(ainx_start + appended_length)
	inx_start = int(numpy.round((mid_freq - grand_spectra_group['axion_frequencies'][0])/(95.4)) - 256/2)
	inx_end = int(inx_start + 256)
	
	if not corrupt_scan and object.op!='nil':
		try:
			#Perform coaddition using class methods
			object.chunk['optimal_weight_sum'][ainx_start:ainx_end] = object.optimal_weight_sum_consolidation(ainx_start, ainx_end)
			object.chunk['model_excess_sqrd'][ainx_start:ainx_end] = object.model_excess_sqrd_consolidation(ainx_start, ainx_end)
			object.chunk['nscans'][ainx_start:ainx_end] = object.nscan_consolidation(ainx_start, ainx_end)
			object.chunk['SNR'][ainx_start:ainx_end] = object.SNR_consolidation(ainx_start, ainx_end)
			object.chunk['axion_fit_uncertainty'][ainx_start:ainx_end] = object.sigma_A_consolidation(ainx_start, ainx_end)
			object.chunk['power_deviation'][inx_start:inx_end] = object.power_deviation_consolidation(inx_start, inx_end) #formerly weighted deltas
			object.chunk['sensitivity_coupling'][ainx_start:ainx_end] = object.coupling_sensitivity_consolidation(ainx_start, ainx_end)
			object.chunk['sensitivity_power'][ainx_start:ainx_end] = object.power_sensitivity_consolidation(ainx_start, ainx_end)
			object.chunk['noise_power'][ainx_start:ainx_end] = object.noise_power_consolidation(ainx_start, ainx_end)
			object.chunk['axion_fit'][ainx_start:ainx_end] = object.axion_fit_consolidation(ainx_start, ainx_end)
			object.chunk['axion_fit_significance'][ainx_start:ainx_end] = object.axion_fit_significance_consolidation(ainx_start, ainx_end)
			object.chunk['deltas'][str(scan_id)][...] = object.scan['deltas']
		except (MemoryError, KeyError, IndexError) as error:
			open(os.getcwd() + '/oo_analysis/meta/error_log', 'a+').write("\n\n"+ str(error))
			print("Error with scan {0} in coaddition script. Writing to error log. ".format(scan_id))
			raise
	
	
	#Update grand spectra last calculation flag
	lastcalc = dt.datetime.now()
	lastcalc = lastcalc.strftime('%Y-%m-%d %H:%M:%S')
	object.chunk['last_change'][0] = str(lastcalc).encode()	
	
	


#Jit is a decorator that saves the function into machine code, speeding up execution
@jit
def divider(self, a1, a2):
	afit = list(range(len(a1)))

	for i in range(len(a1)):
		if a1[i]==0 and a2[i]==0:
			afit[i]=0
		elif a1[i]!=0 and a2[i]==0:
			afit[i]=np.inf
		else:
			afit[i]=a1[i]/a2[i]
	return afit
	
@jit
def quadrature(self, a1, a2, op, pow):
	for i in range(len(a1)):
		if self.op=="+":
			a1[i] = ((a1[i]**pow) + (a2[i]**pow))**(1/2)
		else:
			if (1/a1[i])<(1/a2[i]):
				a1[i]==0
			else:
				a1[i] = ((a1[i]**pow) - (a2[i]**pow))**(1/2)
	return a1

@jit
def inverse_quadrature(self, a1, a2, op, pow):
	for i in range(len(a1)):
		if self.op=="+":
			a1[i] = ((a1[i]**pow) + (a2[i]**pow))**(-1/2)
		else:
			if (1/a1[i])<(1/a2[i]):
				a1[i]==0
			else:
				a1[i] = ((a1[i]**pow) - (a2[i]**pow))**(-1/2)
	return a1

@jit
def inverse_root_quadrature(self, a1, a2, op, pow):
	for i in range(len(a1)):
		if self.op=="+":
			a1[i] = ((a1[i]**pow) + (a2[i]**pow))**(-1/4)
		else:
			if (1/a1[i])<(1/a2[i]):
				a1[i]==0
			else:
				a1[i] = ((a1[i]**pow) - (a2[i]**pow))**(-1/4)
	return a1

@jit
def add_sub(self, a1, a2, op):
	for i in range(len(a1)):
		if op=="+":
			a1[i]+= a2[i]
		else:
			a1[i]-+a2[i]
	return a1






#################### Trash Yard #######################
"""
			
		
			
	
	def frequency_index_matcher(self):
"""
		#Function returns the indices of the scan and chunk where the frequencies match, up to a tolerance. First index is chunk frequencies and second index is scan frequencies
"""
		tol =  49# Hz
		cfreqs = self.chunk_axion_frequencies
		sfreqs = self.scan['axion_frequencies']
		two_indices = np.array([[],[]]) #first array is indices in chunk, second array is indices in scan
		
		si=None
		ci=None
		try:
			if len(cfreqs)==0:
				#if cfreqs is empty (initialized), return zero matching indices
				return two_indices
			
			cmin = np.abs(sfreqs-cfreqs[0])
			if cmin.min()<=tol:		#if sfreqs intersects beginning of cfreqs, set si equal to index where beginnin of cfreq intersect sfreqs
				si = cmin.argmin()
				
			if si==None:
					#if cfreqs intersects beginning of sfreqs, set ci equal to index where beginning of scan intersects cfreqs
				smin = np.abs(cfreqs-sfreqs[0])
				if smin.min()<=tol:
					ci = smin.argmin()
				
						
			if si!=None:
				for sidx in range(len(sfreqs)):
					#Iterate over sfreqs indices and insert the pair of indices, whose frequencies are within 'tol' of each other, into matching_indices
					if (sidx-si)<0:
						pass #Do nothing, since out of range of cfreq
					elif not ((isinstance(sfreqs[sidx], float) or isinstance(sfreqs[sidx], int)) or  (isinstance(cfreqs[sidx-si+1], float) or isinstance(cfreqs[sidx-si+1], int))):
						#do nothing
						pass
					else:
						try:
							if np.abs(sfreqs[sidx]-cfreqs[sidx-si])<=tol:
								two_indices = np.insert(two_indices, 0, [sidx-si, sidx], axis=1)
						except (IndexError, KeyError):
							pass
			elif ci!=None:
				for sidx in np.arange(len(sfreqs)):
					if (len(cfreqs)-1-ci+sidx)<0:
						pass #Do nothing, since out of range of cfreq
					elif not ((isinstance(sfreqs[sidx], float) or isinstance(sfreqs[sidx], int)) or (isinstance(cfreqs[sidx+ci], float) or isinstance(cfreqs[sidx+ci], int))):
						#do nothing
						pass
					elif np.abs(sfreqs[sidx] - cfreqs[sidx+ci])<=tol:
						two_indices = np.insert(two_indices, 0, [sidx+ci, sidx], axis=1)
			elif ci==None and si==None:
				return two_indices
		except (IndexError, KeyError) as error:
			print('frequency matcher failed at scan {0} with error: {1}'.format(self.scan_scan_number, error))
			open(os.getcwd() + '/oo_analysis/meta/error_log', 'a+').write("\n\n"+ str(error))
			raise

		
		#reorder nested lists into accending order
		two_indices[0].sort(); two_indices[1].sort()
		#first index is cfreqs, second is sfreqs
		two_indices = [list(map(int, two_indices[0])), list(map(int, two_indices[1]))] # convert elements to integers for indexing
		return numpy.asarray(two_indices)
		
		
		
		
		
		
def initialize_datapoints(chunk, scan, matched_indices):

	#Determine the frequencies and indices of the scan that dont intersect the chunk frequencies.
	sfreqs = scan["axion_frequencies"]
	matched_freqs = [sfreqs[i] for i in matched_indices[1]]
	if len(matched_freqs)>0:
		aligned_freqs = matched_freqs
	else:
		aligned_freqs=[]
		
	naligned_freqs = []
	leftnaligned_indices = []
	rightnaligned_indices = []
	leftBreak=False
	
	#The following loop isnt very elegant but I need to determine on what side of the chunk does the scan hang off of in order to properly size the chunk.
	for i in sfreqs:
		if i not in aligned_freqs:
			naligned_freqs.append(i)
			if leftBreak==False:
				leftnaligned_indices.append(sfreqs.index(i))	
			elif leftBreak==True:
				rightnaligned_indices.append(sfreqs.index(i))
				
		if i in aligned_freqs:
			leftBreak=True	
	left_hanging = False
	right_hanging = False
	if 0 in leftnaligned_indices:
		#Does the scan frequencies hang off the left side of the chunk frequencies
		left_hanging = True
	if (len(sfreqs)-1) in rightnaligned_indices:
		#Does the scan frequencies hang off the right side of the chunk frequencies
		right_hanging = True
	if left_hanging:
		pos = 0
		for i in leftnaligned_indices[::-1]:
			addtodataset(chunk["sensitivity_coupling"], np.inf, position=pos)
			addtodataset(chunk["axion_fit"], 0, position=pos)
			addtodataset(chunk["axion_fit_uncertainty"], np.inf, position=pos)
			addtodataset(chunk['model_excess_sqrd'], 0, position=pos)
			addtodataset(chunk['nscans'], 0, position=pos)
			addtodataset(chunk['optimal_weight_sum'], 0, position=pos)
			addtodataset(chunk['power_deviation'], 0, position=pos)
			addtodataset(chunk['sensitivity_power'], np.inf, position=pos)
			addtodataset(chunk['SNR'], 0, position=pos)
			addtodataset(chunk['axion_frequencies'], sfreqs[i], position=pos)
	if right_hanging:
		pos=len(sfreqs)-1
		for i in rightnaligned_indices:
			addtodataset(chunk["sensitivity_coupling"], np.inf, position=pos)
			addtodataset(chunk["axion_fit"], 0, position=pos)
			addtodataset(chunk["axion_fit_uncertainty"], np.inf, position=pos)
			addtodataset(chunk['model_excess_sqrd'], 0, position=pos)
			addtodataset(chunk['nscans'], 0, position=pos)
			addtodataset(chunk['optimal_weight_sum'], 0, position=pos)
			addtodataset(chunk['power_deviation'], 0, position=pos)
			addtodataset(chunk['sensitivity_power'], np.inf, position=pos)
			addtodataset(chunk['SNR'], 0, position=pos)
			addtodataset(chunk['axion_frequencies'], sfreqs[i], position=pos)
		
	return chunk
"""