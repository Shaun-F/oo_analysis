
import datetime as dt
import numpy as np
import sys; sys.path.append("..")
from toolbox.add_to_dataset import * #addtodataset, subtractfromdataset, assign_newdata
#import numba; from numba import njit, jit, jitclass, int32, float32
import time
from toolbox.grand_spectra_init import initialization


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
		
		init_scans_in_class_start = time.time()
		self.op = op
		self.scan_scan_number = scid
		self.scan_nscans = scan["nscans"]
		self.scan_sigma_w = scan["sigma_w"]
		self.scan_optimal_weight_sum = scan["optimal_weight_sum"]
		self.scan_SNR = scan["SNR"]
		self.scan_noise_power = scan["noise_power"]
		self.scan_power_deviation = scan["power_deviation"]
		self.scan_model_excess_sqrd = scan["model_excess_sqrd"]
		self.scan_axion_fit = scan["axion_fit"]
		self.scan_axion_fit_uncertainty = scan["axion_fit_uncertainty"]
		self.scan_sensitivity_power = scan["sensitivity_power"]
		self.scan_sensitivity_coupling = scan["sensitivity_coupling"]
		self.scan_axion_frequencies = scan["axion_frequencies"]
		init_scans_in_class_stop = time.time()
		#setup class definitions for chunk 
		
		init_GS_start = time.time()
		#initialization(chunk) #Check if grand spectra group needs to be created
		init_GS_stop = time.time()
		
		#pull only array of points from grand_spectra that line up with scan
		init_GS_in_class_start = time.time()
		self.chunk_axion_frequencies = chunk["axion_frequencies"][...]
		self.matched_indices = self.frequency_index_matcher()
		self.chunk_indices = list(self.matched_indices[0])
		min_inx = min(self.chunk_indices)
		max_inx = max(self.chunk_indices)+1
		self.chunk_scan_number = chunk["scans"][...]
		self.chunk_nscans = chunk["nscans"][min_inx:max_inx]
		self.chunk_sigma_w = chunk["sigma_w"][min_inx:max_inx]
		self.chunk_optimal_weight_sum = chunk["optimal_weight_sum"][min_inx:max_inx]
		self.chunk_SNR = chunk["SNR"][min_inx:max_inx]
		self.chunk_noise_power = chunk["noise_power"][min_inx:max_inx]
		self.chunk_power_deviation = chunk["power_deviation"][min_inx:max_inx]
		self.chunk_model_excess_sqrd = chunk["model_excess_sqrd"][min_inx:max_inx]
		self.chunk_axion_fit = chunk["axion_fit"][min_inx:max_inx]
		self.chunk_axion_fit_uncertainty = chunk["axion_fit_uncertainty"][min_inx:max_inx]
		self.chunk_sensitivity_power = chunk["sensitivity_power"][min_inx:max_inx]
		self.chunk_sensitivity_coupling = chunk["sensitivity_coupling"][min_inx:max_inx]
		self.chunk_axion_fit_significance = chunk['axion_fit_significance'][min_inx:max_inx]
		
		init_GS_in_class_stop = time.time()
		
		submeta = params['submeta']
		if submeta['timeit']:
			submeta['init_scans_in_class'].append(init_scans_in_class_stop - init_scans_in_class_start)
			submeta['init_grand_spectra'].append(init_GS_stop - init_GS_start)
			submeta['init_grand_spectra_in_class'].append(init_GS_in_class_stop-init_GS_in_class_start)
			"""
			if self.op=="+":
				submeta['growing_grand_spectra'].append(growing_GS_stop-growing_GS_start)
				submeta['reinit_grand_spectra_in_class'].append(reinit_GS_in_class_stop-reinit_GS_in_class_start)
			"""
			
			
			
	
	def frequency_index_matcher(self):
		"""
		Function returns the indices of the scan and chunk where the frequencies match, up to a tolerance. First index is chunk frequencies and second index is scan frequencies
		"""
		tol =  49# Hz
		cfreqs = self.chunk_axion_frequencies
		sfreqs = self.scan_axion_frequencies
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
			raise
		
		#reorder nested lists into accending order
		two_indices[0].sort(); two_indices[1].sort()
		#first index is cfreqs, second is sfreqs
		two_indices = [list(map(int, two_indices[0])), list(map(int, two_indices[1]))] # convert elements to integers for indexing
		return numpy.asarray(two_indices)
	
	def axion_fit_consolidation(self):
		optimal_weight_sum = self.chunk_optimal_weight_sum
		model_excess_sqrd = self.chunk_model_excess_sqrd
		axion_fit=list(range(len(optimal_weight_sum)))
		for i in range(len(optimal_weight_sum)):
			if optimal_weight_sum[i]==0 and model_excess_sqrd[i]==0:
				axion_fit[i]=0
			elif optimal_weight_sum[i]!=0 and model_excess_sqrd[i]==0:
				axion_fit[i]=np.inf
			else:
				axion_fit[i]=optimal_weight_sum[i]/model_excess_sqrd[i]
		return axion_fit
		
	def axion_fit_significance_consolidation(self):
		AF = self.chunk_axion_fit
		sigma_A = self.chunk_axion_fit_uncertainty
		AFS = list(range(len(AF)))
		for i in range(len(AFS)):
			if AF[i]==0 and sigma_A[i]==0:
				AFS[i]==0
			elif AF[i]!=0 and sigma_A[i]==0:
				AFS[i]=np.inf
			else:
				AFS[i] = AF[i]/sigma_A[i]
		return AFS
	
	def coupling_sensitivity_consolidation(self):
		csensitivity = self.chunk_sensitivity_coupling
		ssensitivity = self.scan_sensitivity_coupling
		for i, val in enumerate(csensitivity):
			if self.op =="+":
				csensitivity[i] = (1/val**4 + 1/ssensitivity[i]**4)**(-0.25)
			else:
				if (1/val)<(1/ssensitivity[i]):
					csensitivity[i]=0
				else:
					csensitivity[i] = (1/val**4 - 1/ssensitivity[i]**4)**(-0.25)
		return csensitivity
	
	def maximum_likelihood_uncertainty_consolidation(self):
		axion_fit_uncertainty = self.chunk_axion_fit_uncertainty
		scan_axion_fit_uncertainty = self.scan_axion_fit_uncertainty
		for i, val in enumerate(axion_fit_uncertainty):
			if self.op == "+":
				axion_fit_uncertainty[i] = (1/val**2 + 1/scan_axion_fit_uncertainty[i]**2)**(0.5)
			else:
				if (1/val)<(1/scan_axion_fit_uncertainty[i]):
					axion_fit_uncertainty[i]=0
				else:
					axion_fit_uncertainty[i] = (1/val**2 - 1/scan_axion_fit_uncertainty[i]**2)**(0.5)
		return axion_fit_uncertainty
	
	def model_excess_sqrd_consolidation(self):
		MES = self.chunk_model_excess_sqrd
		sMES = self.scan_model_excess_sqrd
		for i, val in enumerate(MES):

			if self.op=="+":
				MES[i] = val + sMES[i]
			else:
				MES[i] = val - sMES[i]
		return MES
	
	def noise_power_consolidation(self):
		SNR = self.chunk_SNR
		WS = self.chunk_optimal_weight_sum
		noise_power=np.asarray([])
		if len(np.where(WS==0))!=0:
			for inx, val in enumerate(WS):
				if val==0 and SNR[inx]==0:
					#THis could be the case when grand spectra bin is initialized. Maybe theres a better way of doing this?
					noise_power = np.append(noise_power, 0)
				elif val==0 and SNR[inx]!=0:
					noise_power = np.append(noise_power, np.inf)
				else:
					noise_power = np.append(noise_power, SNR[inx]/val)
		else:
			noise_power=SNR/WS
		return noise_power
	
	def nscan_consolidation(self):
		nscans=self.chunk_nscans
		scan_nscan=self.scan_nscans
		for i, val in enumerate(nscans):
			
			if self.op=="+":
				nscans[i] = val + scan_nscan[i]
			else:
				nscans[i] = val - scan_nscan[i]
		return nscans
		
	def optimal_weight_sum_consolidation(self):
		WS = self.chunk_optimal_weight_sum
		sWS = self.scan_optimal_weight_sum
		for i, val in enumerate(WS):
			if self.op=="+":
				WS[i] = val + sWS[i]
			else:
				WS[i] = val - sWS[i]
		return WS
	
	def power_deviation_consolidation(self):
		power_deviation = self.chunk_power_deviation
		spower_deviation = self.scan_power_deviation
		for i, val in enumerate(power_deviation):
			if self.op=="+":
				power_deviation[i] = val + spower_deviation[i]
			else:
				power_deviation[i] = val - spower_deviation[i]
		return power_deviation
	
	def power_sensitivity_consolidation(self):
		sensitivity = self.chunk_sensitivity_power
		ssensitivity = self.scan_sensitivity_power
		for i, val in enumerate(sensitivity):
			if self.op=="+":
				sensitivity[i]=(1/val**2 + 1/ssensitivity[i]**2)**(-0.5)
			else:
				if (1/val)<(1/ssensitivity[i]):
					sensitivity[i]=0
				else:
					sensitivity[i]=(1/val**2 - 1/ssensitivity[i]**2)**(-0.5)
		return sensitivity
		
	def sigma_A_consolidation(self):
		sigma_A = self.chunk_axion_fit_uncertainty
		ssigma_A = self.scan_axion_fit_uncertainty
		for i, val in enumerate(sigma_A):
			if self.op =="+":
				sigma_A[i] = (1/(val**2) + 1/(ssigma_A[i]**2))**(1/2)
			else:
				if (1/val)<(1/ssigma_A[i]):
					sigma_A[i]=0
				else:
					sigma_A[i] = (np.abs(1/(val**2) - 1/(ssigma_A[i]**2)))**(1/2)
		return sigma_A
	
	def SNR_consolidation(self):
		cSNR = self.chunk_SNR
		sSNR = self.scan_SNR
		for i, val in enumerate(cSNR):
			if self.op == "+":
				cSNR[i] = (val**2 + sSNR[i]**2)**(1/2)
			else:
				if val<sSNR[i]:
					cSNR[i]=0
				else:
					cSNR[i] = (val**2 - sSNR[i]**2)**(1/2)
		return cSNR
	
	def weighted_delta_consolidation(self):
		weighted_deltas = self.chunk_weighted_deltas
		swds = self.scan_weighted_deltas
		for i, val in enumerate(weighted_deltas):

			if self.op=="+":
				weighted_deltas[i] = val + swds[i]
			else:
				weighted_deltas[i] = val - swds[i]
		return weighted_deltas
	
	def close_out(self, chunk):
		chunk['nscans'][self.chunk_indices] = self.chunk_nscans
		chunk['sigma_w'][self.chunk_indices] = self.chunk_sigma_w
		chunk['optimal_weight_sum'][self.chunk_indices] = self.chunk_optimal_weight_sum
		chunk['SNR'][self.chunk_indices] = self.chunk_SNR
		chunk['noise_power'][self.chunk_indices] = self.chunk_noise_power
		chunk['power_deviation'][self.chunk_indices] = self.chunk_power_deviation
		chunk['model_excess_sqrd'][self.chunk_indices] = self.chunk_model_excess_sqrd
		chunk['axion_fit'][self.chunk_indices] = self.chunk_axion_fit
		chunk['axion_fit_uncertainty'][self.chunk_indices] = self.chunk_axion_fit_uncertainty
		chunk['sensitivity_power'][self.chunk_indices] = self.chunk_sensitivity_power
		chunk['sensitivity_coupling'][self.chunk_indices] = self.chunk_sensitivity_coupling
		chunk['scans_in'][...] = chunk['scans_in'][...]
		chunk['scans_out'][...] = chunk['scans_out'][...]
		

		

def add_subtract_scan(add_subtract, scan, chunk, scan_id, **params):
	"""
	Parameters
		add_subtract: ('add', 'subtract') Determines what to do with scan
		scan: (dictionary) items to add to grand spectra, including deltas
		chunk: (h5py dataset) dataset to add to.
		scan_id: (
	"""
	#initialize scan object	
	if add_subtract == "add":
		op = "+"
	elif add_subtract == "subtract":
		op = "-"
	else:
		return "Error: combining operation not recognized. Available operations are 'add' and 'subtract'"
		
	#Running initiliazation procedure for chunks, i.e. create missing arrays.
	#chunk = initialization(chunk) 									Not currently using.
	if op=="+" and (scan_id.encode() in chunk["scans_in"]):
		phrase = "scan " + str(scan_id) + " already added. Skipping coaddition"
		#print(phrase)
		return phrase
	elif op=="-" and (scan_id.encode() in chunk["scans_out"]):
		phrase = "scan " + str(scan_id) + " already subtracted. Skipping coaddition"
		#print(phrase)
		return phrase
	
	
	class_init_start = time.time()
	scan_obj = scan_cl(scan=scan, scid = scan_id, chunk=chunk, op=op, **params)
	class_init_stop = time.time()
	
	corrupt_scan=False
	if type(scan) is str:
		corrupt_scan = True
		add_subtract = 'ommit'
		op = 'nil'


	consolidation_time_start = time.time()
	if not corrupt_scan and op!='nil':
		opt_wght_sum_start = time.time()
		scan_obj.chunk_optimal_weight_sum = scan_obj.optimal_weight_sum_consolidation()
		opt_wght_sum_stop = time.time()
		scan_obj.chunk_model_excess_sqrd = scan_obj.model_excess_sqrd_consolidation()
		nscans_stop = time.time()
		scan_obj.chunk_nscans = scan_obj.nscan_consolidation()
		SNR_stop = time.time()
		scan_obj.chunk_SNR = scan_obj.SNR_consolidation()
		sigma_A_stop = time.time()
		scan_obj.chunk_axion_fit_uncertainty = scan_obj.sigma_A_consolidation()
		power_dev_stop = time.time()
		scan_obj.chunk_power_deviation = scan_obj.power_deviation_consolidation() #formerly weighted deltas
		coupl_sens_stop = time.time()
		scan_obj.chunk_sensitivity_coupling = scan_obj.coupling_sensitivity_consolidation()
		power_sens_stop = time.time()
		scan_obj.chunk_sensitivity_power = scan_obj.power_sensitivity_consolidation()
		noise_pow_stop = time.time()
		scan_obj.chunk_noise_power = scan_obj.noise_power_consolidation()
		axion_fit_stop = time.time()
		scan_obj.chunk_axion_fit = scan_obj.axion_fit_consolidation()
		axion_fit_sig_stop = time.time()
		scan_obj.chunk_axion_fit_significance = scan_obj.axion_fit_significance_consolidation()
		scans_start = time.time()
		addtodataset(chunk['scans'], scan_id)
		scans_stop = time.time()
	consolidation_time_stop = time.time()
	
	scans_in_start = time.time()
	if add_subtract=='add':
		addtodataset(chunk['scans_in'], scan_id)
		subtractfromdataset(chunk['scans_out'], array_or_string=scan_id)
	elif add_subtract=="subtract":
		addtodataset(chunk['scans_out'], scan_id)
		subtractfromdataset(chunk['scans_in'], array_or_string=scan_id)
	scans_in_stop = time.time()
	
	
	lastcalc = dt.datetime.now()
	lastcalc = lastcalc.strftime('%Y-%m-%d %H:%M:%S')
	chunk.attrs["last_change"] = str(lastcalc)

	
	submeta = params['submeta']
	if submeta['timeit']:
		submeta['consolidation_time'].append(consolidation_time_stop-consolidation_time_start)
		submeta['scans_in_time'].append(scans_in_stop-scans_in_start)
		submeta['class_init_time'].append(class_init_stop-class_init_start)
		submeta['optimal_weight_sum_consolidation_time'].append(opt_wght_sum_stop-opt_wght_sum_start)
		submeta['model_excess_sqrd_consolidation_time'].append(nscans_stop-opt_wght_sum_stop)
		submeta['nscan_consolidation_time'].append(SNR_stop-nscans_stop)
		submeta['SNR_consolidation_time'].append(sigma_A_stop-SNR_stop)
		submeta['sigma_A_consolidation_time'].append(power_dev_stop-sigma_A_stop)
		submeta['power_deviation_consolidation_time'].append(coupl_sens_stop-power_dev_stop)
		submeta['coupling_sensitivity_consolidation_time'].append(power_sens_stop-coupl_sens_stop)
		submeta['power_sensitivity_consolidation_time'].append(noise_pow_stop-power_sens_stop)
		submeta['noise_power_consolidation_time'].append(axion_fit_stop-noise_pow_stop)
		submeta['axion_fit_consolidation_time'].append(axion_fit_sig_stop-axion_fit_stop)
		submeta['axion_fit_significance_consolidation_time'].append(scans_start-axion_fit_sig_stop)
		submeta['scans_in_grand_spectra_addition_time'].append(scans_stop-scans_start)
	
	scan_obj.close_out(chunk)
	
	
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
#initialize attribute arrays to be populated (coadded)
if self.op=="+":
	growing_GS_start = time.time()
	chunk = initialize_datapoints(chunk, scan, self.matched_indices) 
	growing_GS_stop = time.time()
	#Initialize any attributes not previously defined
	reinit_GS_in_class_start = time.time()
	self.chunk_axion_frequencies = chunk["axion_frequencies"][...]
	self.matched_indices = self.frequency_index_matcher()
	self.chunk_scan_number = chunk["scans"][...]
	self.chunk_nscans = chunk["nscans"][...]
	self.chunk_sigma_w = chunk["sigma_w"][...]
	self.chunk_optimal_weight_sum = chunk["optimal_weight_sum"][...]
	self.chunk_SNR = chunk["SNR"][...]
	self.chunk_noise_power = chunk["noise_power"][...]
	self.chunk_power_deviation = chunk["power_deviation"][...]
	self.chunk_model_excess_sqrd = chunk["model_excess_sqrd"][...]
	self.chunk_axion_fit = chunk["axion_fit"][...]
	self.chunk_axion_fit_uncertainty = chunk["axion_fit_uncertainty"][...]
	self.chunk_sensitivity_power = chunk["sensitivity_power"][...]
	self.chunk_sensitivity_coupling = chunk["sensitivity_coupling"][...]
	reinit_GS_in_class_stop = time.time()
	
elif self.op=="-":
	print('subtracting')
	subtractfromdataset(chunk['axion_frequencies'], scan['axion_frequencies'])
"""