
import datetime as dt
import numpy as np

class scan_cl(object):
	"""
	This class performs the necessary addition or subtraction of scans from a grand spectrum.
	Attributes: scan_number, ch, res, modulation_type, axion_shape
	Methods: axion_scan, digitizer_scan, analyze_scan, add_scan
	"""
	
	
	def __init__(self, scid, scan, chunk, op):
		"""
		Initialization method
		"""
		self.op = op
		
		self.scan_scan_number = scid
		self.scan_nscans = scan["nscans"]
		self.scan_sigma_w = scan["sigma_w"]
		self.scan_optimal_weight_sum = scan["optimal_weight_sum"]
		self.scan_SNR = scan["SNR"]
		self.scan_noise_power = scan["noise_power"]
		self.scan_weighted_deltas = scan["power_deviation"]
		self.scan_model_excess_sqrd = scan["model_excess_sqrd"]
		self.scan_axion_fit = scan["axion_fit"]
		self.scan_axion_fit_uncertainty = scan["axion_fit_uncertainty"]
		self.scan_sensitivity_power = scan["sensitivity_power"]
		self.scan_sensitivity_coupling = scan["sensitivity_coupling"]
		self.scan_axion_frequencies = scan["axion_frequencies"]
		
		#setup class definitions for chunk 
		initialization(chunk,scan)
		self.chunk_scan_number = chunk.attrs["scans"]
		self.chunk_nscans = chunk.attrs["nscans"]
		self.chunk_sigma_w = chunk.attrs["sigma_w"]
		self.chunk_optimal_weight_sum = chunk.attrs["optimal_weight_sum"]
		self.chunk_SNR = chunk.attrs["SNR"]
		self.chunk_noise_power = chunk.attrs["noise_power"]
		self.chunk_weighted_deltas = chunk.attrs["power_deviation"]
		self.chunk_model_excess_sqrd = chunk.attrs["model_excess_sqrd"]
		self.chunk_axion_fit = chunk.attrs["axion_fit"]
		self.chunk_axion_fit_uncertainty = chunk.attrs["axion_fit_uncertainty"]
		self.chunk_sensitivity_power = chunk.attrs["sensitivity_power"]
		self.chunk_sensitivity_coupling = chunk.attrs["sensitivity_coupling"]
		self.chunk_axion_frequencies = chunk.attrs["axion_frequencies"]
		
		#initialize attribute arrays to be populated (coadded)
		chunk = initialize_datapoints(chunk, scan, self.frequency_index_matcher()) #Initialize any attributes not previously defined
		self.chunk_scan_number = np.sort(chunk.attrs["scans"])
		self.chunk_nscans = np.sort(chunk.attrs["nscans"])
		self.chunk_sigma_w = np.sort(chunk.attrs["sigma_w"])
		self.chunk_optimal_weight_sum = np.sort(chunk.attrs["optimal_weight_sum"])
		self.chunk_SNR = np.sort(chunk.attrs["SNR"])
		self.chunk_noise_power = np.sort(chunk.attrs["noise_power"])
		self.chunk_weighted_deltas = np.sort(chunk.attrs["power_deviation"])
		self.chunk_model_excess_sqrd = np.sort(chunk.attrs["model_excess_sqrd"])
		self.chunk_axion_fit = np.sort(chunk.attrs["axion_fit"])
		self.chunk_axion_fit_uncertainty = np.sort(chunk.attrs["axion_fit_uncertainty"])
		self.chunk_sensitivity_power = np.sort(chunk.attrs["sensitivity_power"])
		self.chunk_sensitivity_coupling = np.sort(chunk.attrs["sensitivity_coupling"])
		self.chunk_axion_frequencies = np.sort(chunk.attrs["axion_frequencies"])
	
	def frequency_index_matcher(self):
		"""
		Function returns the indices of the scan and chunk where the frequencies match, up to a tolerance. First index is chunk frequencies and second index is scan frequencies
		"""
		tol = 1 # Hz
		cfreqs = self.chunk_axion_frequencies
		sfreqs = self.scan_axion_frequencies
		
		si=None
		ci=None
		if len(cfreqs)==0:
			#if cfreqs is empty (initialized), return zero matching indices
			return []
			
		for inx, val in enumerate(sfreqs):
			if np.abs(cfreqs[0]-val)<=tol:
				#if sfreqs intersects beginning of cfreqs, set si equal to index where beginnin of cfreq intersect sfreqs
				si = inx
				break
		if si==None:
				#if cfreqs intersects beginning of sfreqs, set ci equal to index where beginning of scan intersects cfreqs
			for cidx, cfreq in enumerate(cfreqs):
				if np.abs(sfreqs[0]-cfreq)<=tol:
					ci = cidx
					break
					
		two_indices = np.array([[],[]])
		if si!=None:
			for sidx in range(len(sfreqs)):
				#Iterate over sfreqs indices and insert the pair of indices, whose frequencies are within 'tol' of each other, into matching_indices
				if (sidx-si)<0:
					pass #Do nothing, since out of range of cfreq
				elif not (isinstance(sfreqs[sidx], float) and isinstance(sfreqs[sidx], int)) or not (isinstance(cfreqs[sidx-si+1], float) and isinstance(cfreqs[sidx-si+1], int)):
					#do nothing
					pass
				elif float(sfreqs[sidx]-cfreqs[sidx-si])<=tol:
					np.insert(two_indices, 0, [sidx-si,sidx], axis=1)
		elif ci!=None:
			for sidx in np.arange(len(sfreqs)):
				if (np.argwhere(cfreqs==cfreqs[-1])[0][0]-ci-sidx)<0:
					pass #Do nothing, since out of range of cfreq
				elif not (isinstance(sfreqs[sidx-ci+1], float) and isinstance(sfreqs[sidx-ci+1], int)) or not (isinstance(cfreqs[sidx], float) and isinstance(cfreqs[sidx], int)):
					#do nothing
					pass
				elif np.abs(sfreqs[sidx] - cfreqs[sidx+ci])<=tol:
					np.insert(two_indices, 0, [sidx+ci, sidx], axis=1)
		else:
			return two_indices
		#reorder nested lists into accending order
		two_indices[0].sort(); two_indices[1].sort()
		
		#first index is cfreqs, second is sfreqs
		return two_indices
	
	def axion_fit_consolidation(self):
		optimal_weight_sum = self.chunk_optimal_weight_sum
		model_excess_sqrd = self.chunk_model_excess_sqrd
		axion_fit=list(range(len(optimal_weight_sum)))
		for i in range(len(optimal_weight_sum)):
			if optimal_weight_sum[i]==0 and model_excess_sqrd[i]==0:
				axion_fit[i]=0
			elif optimal_weight_sum[i]!=0 and model_excess_sqrd[i]==0:
				axion_fit[i]=float("inf")
			else:
				axion_fit[i]=optimal_weight_sum[i]/model_excess_sqrd[i]
		return axion_fit
		
	def axion_fit_significance_consolidation(self):
		AF = self.chunk_axion_fit
		sigma_A = self.chunk_axion_fit_uncertainty
		AFS = AF/sigma_A
		return AFS
	
	def coupling_sensitivity_consolidation(self):
		csensitivity = self.chunk_sensitivity_coupling
		ssensitivity = self.scan_sensitivity_coupling
		matching_indices = self.frequency_index_matcher()
		cidx = matching_indices[0]
		sidx = matching_indices[1]
		for i, val in enumerate(cidx):
			if self.op =="+":
				csensitivity[val] = (1/csensitivity[val]**4 + 1/ssensitivity[sidx[i]]**4)**(-0.25)
			else:
				csensitivity[val] = (1/csensitivity[val]**4 - 1/ssensitivity[sidx[i]]**4)**(-0.25)
		return csensitivity
	
	def maximum_likelihood_uncertainty_consolidation(self):
		matching_indices = self.frequency_index_matcher()
		axion_fit_uncertainty = self.chunk_axion_fit_uncertainty
		scan_axion_fit_uncertainty = self.scan_axion_fit_uncertainty
		cidx = matching_indices[0]
		sidx = matching_indices[1]
		for i, val in enumerate(cidx):
			if self.op == "+":
				axion_fit_uncertainty[val] = (1/axion_fit_uncertainty[val]**2 + 1/scan_axion_fit_uncertainty[sidx[i]]**2)**(0.5)
			else:
				axion_fit_uncertainty[val] = (1/axion_fit_uncertainty[val]**2 - 1/scan_axion_fit_uncertainty[sidx[i]]**2)**(0.5)
		return axion_fit_uncertainty
	
	def model_excess_sqrd_consolidation(self):
		matching_index = self.frequency_index_matcher()
		MES = self.chunk_model_excess_sqrd
		sMES = self.scan_model_excess_sqrd
		cidx = matching_index[0]
		sidx = matching_index[1]
		for i, val in enumerate(cidx):

			if op=="+":
				MES[val] = MES[val] + sMES[sidx[i]]
			else:
				MES[val] = MES[val] - sMES[sidx[i]]
		return MES
	
	def noise_power_consolidation(self):
		SNR = self.chunk_SNR
		WS = self.chunk_optimal_weight_sum
		noise_power=np.asarray([])
		if len(np.where(WS==0))!=0:
			for inx, val in enumerate(WS):
				if val==0:
					#THis could be the case when grand spectra bin is initialized. Maybe theres a better way of doing this?
					noise_power = np.append(noise_power, SNR[inx])
				else:
					noise_power = np.append(noise_power, SNR[inx]/val)
		else:
			noise_power=SNR/WS
		return noise_power
	
	def nscan_consolidation(self):
		indices = self.frequency_index_matcher()
		nscans=self.chunk_nscans
		scan_nscan=self.scan_nscans
		cidx = indices[0]
		sidx = indices[1]
		for i, val in enumerate(cidx):
			
			if self.op=="+":
				nscans[val] = nscans[val] + scan_nscans[sidx[i]]
			else:
				nscans[val] = nscans[val] - scan_nscans[sidx[i]]
		return nscans
		
	def optimal_weight_sum_consolidation(self):
		indices = self.frequency_index_matcher()
		WS = self.chunk_optimal_weight_sum
		sWS = self.scan_optimal_weight_sum
		cidx = indices[0]
		sidx = indices[1]
		for i, val in enumerate(cidx):
			if self.op=="+":
				WS[val] = WS[val] + sWS[sidx[i]]
			else:
				WS[val] = WS[val] - sWS[sidx[i]]
		return WS
	
	def power_deviation_consolidation(self):
		matching_indices = self.frequency_index_matcher()
		power_deviation = self.chunk_power_deviation
		spower_deviation = self.scan_power_deviation
		cidx = indices[0]
		sidx = indices[1]
		for i, val in enumerate(cidx):

			if self.op=="+":
				power_deviation[val] = power_deviation[val] + spower_deviation[sidx[i]]
			else:
				power_deviation[val] = power_deviation[val] - spower_deviation[sidx[i]]
		return power_deviation
	
	def power_senstivity_consolidation(self):
		indices = self.frequency_index_matcher()
		sensitivity = self.chunk_sensitivity_power
		ssensitivity = self.scan_sensitivity_power
		cidx = indices[0]
		sidx = indices[1]
		for i, val in enumerate(cidx):
	
			if self.op=="+":
				sensitivity[val]=(1/sensitivity[val]**2 + 1/ssensitivity[sidx[i]]**2)**(-0.5)
			else:
				sensitivity[val]=(1/sensitivity[val]**2 - 1/ssensitivity[sidx[i]]**2)**(-0.5)
		return sensitivity
		
	def sigma_A_consolidation(self):
		matching_indices = self.frequency_index_matcher()
		sigma_A = self.chunk_axion_fit_uncertainty
		ssigma_A = self.scan_axion_fit_uncertainty
		cidx = matching_indices[0]
		sidx = matching_indices[1]
		for i, val in enumerate(cidx):
			if self.op =="+":
				sigma_A[val] = (1/(sigma_A[val]**2) + 1/(ssigma_A[sidx[i]]**2))**(1/2)
			else:
				sigma_A[val] = (1/(sigma_A[val]**2) - 1/(ssigma_A[sidx[i]]**2))**(1/2)
		return sigma_A
	
	def SNR_consolidation(self):
		matching_indices = self.frequency_index_matcher()
		cSNR = self.chunk_SNR
		sSNR = self.scan_SNR
		cidx = matching_indices[0]
		sidx = matching_indices[1]
		for i, val in enumerate(cidx):

			if self.op == "+":
				cSNR[val] = (cSNR[val]**2 + sSNR[sidx[i]]**2)**(1/2)
			else:
				cSNR[val] = (cSNR[val]**2 - sSNR[sidx[i]]**2)**(1/2)
		return cSNR
	
	def weighted_delta_consolidation(self):
		matching_indices = self.frequency_index_matcher()
		weighted_deltas = self.chunk_weighted_deltas
		swds = self.scan_weighted_deltas
		cidx = matching_indices[0]
		sidx = matching_indices[1]
		for i, val in enumerate(cidx):

			if self.op=="+":
				weighted_deltas[val] = weighted_deltas[val] + swds[sidx[i]]
			else:
				weighted_deltas[val] = weighted_deltas[val] - swds[sidx[i]]
		return weighted_deltas

		
def initialization(chunk, scan):
	#Initialize all attributes with empty containers
	if 'scans' not in chunk.attrs:
		chunk.attrs['scans'] = np.asarray([], dtype='byte') #array of scan id strings
	if 'sigma_w' not in chunk.attrs:
		chunk.attrs['sigma_w'] = []
	if 'optimal_weight_sum' not in chunk.attrs:
		chunk.attrs['optimal_weight_sum']=[]
	if 'SNR' not in chunk.attrs:
		chunk.attrs['SNR'] = []
	if 'noise_power' not in chunk.attrs:
		chunk.attrs['noise_power'] = []
	if 'model_excess_sqrd' not in chunk.attrs:
		chunk.attrs['model_excess_sqrd'] = []
	if 'axion_fit' not in chunk.attrs:
		chunk.attrs['axion_fit'] = []
	if 'axion_fit_uncertainty' not in chunk.attrs:
		chunk.attrs['axion_fit_uncertainty'] = []
	if 'sensitivity_power' not in chunk.attrs:
		chunk.attrs['sensitivity_power'] = []
	if 'sensitivity_coupling' not in chunk.attrs:
		chunk.attrs['sensitivity_coupling'] = []
	if 'axion_frequencies' not in chunk.attrs:
		chunk.attrs['axion_frequencies'] = []
	if 'power_deviation' not in chunk.attrs:
		chunk.attrs["power_deviation"] = []
	if 'nscans' not in chunk.attrs:
		chunk.attrs['nscans'] = []
	if 'scans_in' not in chunk.attrs:
		chunk.attrs['scans_in'] = np.asarray([], dtype='byte') #array of scan id strings in grand spectra
	if 'scans_out' not in chunk.attrs:
		chunk.attrs['scans_out'] = np.asarray([], dtype='byte') #array of scan id strings not in grand spectra
	if 'last_change' not in chunk.attrs:
		chunk.attrs['last_change'] = np.asarray([], dtype='byte') #array of timestamp strings
def initialize_datapoints(chunk, scan, matched_freqs):

	#Determine the frequencies and indices of the scan that dont intersect the chunk frequencies.
	sfreqs = scan["axion_frequencies"]
	
	if len(matched_freqs[1])>0:
		print("Matching indices", matched_freqs)
		aligned_freqs = matched_freqs[1]
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
		for i in leftnaligned_indices:
			chunk.attrs["sensitivity_coupling"]=np.insert(chunk.attrs["sensitivity_coupling"], 0, 0)
			chunk.attrs["axion_fit"]=np.insert(chunk.attrs["axion_fit"], 0, 0)
			chunk.attrs["axion_fit_uncertainty"]=np.insert(chunk.attrs["axion_fit_uncertainty"], 0, np.inf) 
			chunk.attrs["model_excess_sqrd"]=np.insert(chunk.attrs["model_excess_sqrd"], 0, 0)
			chunk.attrs["nscans"]=np.insert(chunk.attrs["nscans"], 0, 0)
			chunk.attrs["optimal_weight_sum"]=np.insert(chunk.attrs["optimal_weight_sum"], 0, 0)
			chunk.attrs["power_deviation"]=np.insert(chunk.attrs["power_deviation"], 0,0)
			chunk.attrs["sensitivity_power"]=np.insert(chunk.attrs["sensitivity_power"], 0, float("inf"))
			chunk.attrs["SNR"]=np.insert(chunk.attrs["SNR"], 0, 0)
			chunk.attrs["axion_frequencies"]=np.insert(chunk.attrs["axion_frequencies"], 0, sfreqs[i])
	if right_hanging:
		for i in rightnaligned_indices:
			chunk.attrs["sensitivity_coupling"]=np.insert(chunk.attrs["sensitivity_coupling"], (len(sfreqs)-1), 0)
			chunk.attrs["axion_fit"]=np.insert(chunk.attrs["axion_fit"], (len(sfreqs)-1), 0)
			chunk.attrs["axion_fit_uncertainty"]=np.insert(chunk.attrs["axion_fit_uncertainty"], (len(sfreqs)-1), float("inf")) 
			chunk.attrs["model_excess_sqrt"]=np.insert(chunk.attrs["model_excess_sqrt"], (len(sfreqs)-1), 0)
			chunk.attrs["nscans"]=np.insert(chunk.attrs["nscans"], (len(sfreqs)-1), 0)
			chunk.attrs["optimal_weight_sum"]=np.insert(chunk.attrs["optimal_weight_sum"], (len(sfreqs)-1), 0)
			chunk.attrs["power_deviation"]=np.insert(chunk.attrs["power_deviation"], (len(sfreqs)-1),0)
			chunk.attrs["sensitivity_power"]=np.insert(chunk.attrs["sensitivity_power"], (len(sfreqs)-1), float("inf"))
			chunk.attrs["SNR"]=np.insert(chunk.attrs["SNR"], (len(sfreqs)-1), 0)
			chunk.attrs["axion_frequencies"]=np.insert(chunk.attrs["axion_frequencies"], (len(sfreqs)-1), sfreqs[i])
		
	return chunk
		

def add_subtract_scan(add_subtract, scan, chunk, scan_id):
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
	scan_obj = scan_cl(scan=scan, scid = scan_id, chunk=chunk, op="+")
	
	#Running initiliazation procedure for chunks, i.e. create missing arrays.
	#chunk = initialization(chunk) 									Not currently using.
	
	if op=="+" and (scan_id in chunk.attrs["scans_in"]):
		phrase = "scan " + str(scan_id) + " already added"
		print(phrase)
		return phrase
	elif op=="-" and (scan_id in chunk.attrs["scans_out"]):
		phrase = "scan " + str(scan_id) + " already subtracted"
		print(phrase)
		return phrase
	
	corrupt_scan=False
	if type(scan) is str:
		corrupt_scan = True
		add_subtract = 'ommit'
		op = nil
	if not corrupt_scan:
		chunk.attrs["nscans"] = scan_obj.nscan_consolidation()
		chunk.attrs["SNR"] = scan_obj.SNR_consolidation()
		chunk.attrs["optimal_weight_sum"] = scan_obj.optimal_weight_sum_consolidation()
		chunk.attrs["noise_power"] = scan_obj.noise_power_consolidation()
		chunk.attrs["model_excess_sqrd"] = scan_obj.model_excess_sqrd_consolidation()
		chunk.attrs["axion_fit_uncertainty"] = scan_obj.sigma_A_consolidation()
		chunk.attrs["axion_fit"] = scan_obj.axion_fit_consolidation()
		chunk.attrs["power_deviation"] = scan_obj.weighted_delta_consolidation() #formerly weighted deltas
		chunk.attrs["axion_fit_significance"] = scan_obj.axion_fit_significance_consolidation()
		chunk.attrs["sensitivity_coupling"] = scan_obj.coupling_sensitivity_consolidation()
		chunk.attrs["sensitivity_power"] = scan_obj.power_senstivity_consolidation()
	chunk.attrs["scans"] = np.append(chunk.attrs["scans"], scan_id).astype('S')

	print(chunk.attrs["power_deviation"])
	if add_subtract=='add':
		chunk.attrs["scans_in"] = np.append(chunk.attrs["scans_in"],scan_id).astype('S')
		for i, chunk_scan_id in enumerate(chunk.attrs["scans_out"]):
			if scan_id == chunk_scan_id:
				chunk.attrs["scans_out"][i]=nil
	elif add_subtract=="subtract":
		chunk.attrs["scans_out"] = np.append(chunk.attrs["scans_out"],scan_id).astype('S')
		for i, chunk_scan_id in enumerate(chunk.attrs["scans_in"]):
			if scan_id == chunk_scan_id:
				chunk.attrs["scans_id"][i] = nil
	
	lastcalc = dt.datetime.now()
	lastcalc = lastcalc.strftime('%Y-%m-%d %H:%M:%S')
	chunk.attrs["last_change"] = str(lastcalc)
	