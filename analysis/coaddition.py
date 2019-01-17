
import datetime as dt
import numpy as np

class scan(object):
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
		self.scan_sigma_w = scan["simga_w"]
		self.scan_optimal_weight_sum = scan["optimal_weigth_sum"]
		self.scan_SNR = scan["SNR"]
		self.scan_noise_power = scan["noise_power"]
		self.scan_weighted_deltas = scan["weighted_deltas"]
		self.scan_model_excess_sqrd = scan["model_excess_sqrd"]
		self.scan_axion_fit = scan["axion_fit"]
		self.scan_axion_fit_uncertainty = scan["axion_fit_uncertainty"]
		self.scan_sensitivity_power = scan["sensitivity_power"]
		self.scan_sensitivity_coupling = scan["sensitivity_coupling"]
		self.scan_axion_frequencies = scan["axion_frequencies"]
		self.scan_power_deviation = scan["power_deviation"]
		
		chunk = initialization(chunk) #Initialize any attributes not previously defined
		
		self.chunk_scan_number = chunk.attrs["scans"]
		self.chunk_nscans = chunk.attrs["nscans"]
		self.chunk_sigma_w = chunk.attrs["sigma_w"]
		self.chunk_optimal_weight_sum = chunk.attrs["optimal_weight_sum"]
		self.chunk_SNR = chunk.attrs["SNR"]
		self.chunk_noise_power = chunk.attrs["noise_power"]
		self.chunk_weighted_deltas = chunk.attrs["weighted_deltas"]
		self.chunk_model_excess_sqrd = chunk.attrs["model_excess_sqrd"]
		self.chunk_axion_fit = chunk.attrs["axion_fit"]
		self.chunk_axion_fit_uncertainty = chunk.attrs["axion_fit_uncertainty"]
		self.chunk_sensitivity_power = chunk.attrs["sensitivity_power"]
		self.chunk_sensitivity_coupling = chunk.attrs["sensitivity_coupling"]
		self.chunk_axion_frequencies = chunk.attrs["axion_frequencies"]
		self.chunk_power_deviation = chunk.attrs["power_deviation"]
	
	
	def frequency_index_matcher(self):
		tol = 0.1
		cfreqs = self.chunk_axion_frequencies
		sfreqs = self.scan_axion_frequencies
		
		for inx, val in enumerate(sfreqs):
			if np.abs(cfreqs[1]-val)<=tol:
				si = sidx
				break
		if not si:
			for cidx, cfreq in enumerate(cfreqs):
				if np.abs(sfreqs[1]-cfreq)<=tol:
					ci = cidx
					break
					
		two_indices = []
		if si:
			for i in np.arange(len(sfreqs)):
				if (type(sfreqs[cidx]) != float and type(sfreqs[cidx])!=int) or (type(cfreqs[idx-si+1])!=float and type(cfreqs[cidx-si+1]))!=int:
					#do nothing
					pass
				elif float(sfreqs[idx]-cfreqs[idx-si+1])<=tol:
					np.insert(two_indices, [idx-si+1,idx])
		elif ci:
			for i in np.arange(len(sfreqs)+ci):
				if type(sfreqs[idx-ci+1])!=float and type(sfreqs[idx-ci+1])!=int or type(cfreqs[idx])!=float and type(cfreqs[idx])!=int:
					#do nothing
					pass
				elif np.abs(sfreqs[idx-ci+1] - cfreqs[idx])<=tol:
					np.insert(two_indices, [idx, idx-ci+1])
		else:
			return "Error: frequencies of scan and chunk do not match"
		return two_indices
	
	def axion_fit_consolidation(self):
		optimal_weight_sum = self.chunk_optimal_weight_sum
		model_excess_sqrd = self.chunk_model_excess_sqrd
		axion_fit = optimal_weight_sum/model_excess_sqrd
		return axion_fit
		
	def axion_fit_significance_consolidation(self):
		AF = self.chunk_axion_fit
		sigma_A = self.chunk_axion_fit_uncertainty
		AFS = AF/sigma_A
		return AFS
	
	def coupling_senstivity_consolidation(self):
		csensitivity = self.chunk_sensitivity_coupling
		ssensitivity = self.scan_sensitivity_coupling
		matching_indices = frequency_index_matcher(self)
		for i, val in enumerate(matching_indices):
			cidx = indices[0]
			sidx = indices[1]
			if self.op =="+":
				sensitivity[cidx] = (1/csensitivity[cidx]**4 + 1/ssensitivity[sidx]**4)**(-0.25)
			else:
				sensitivity[cidx] = (1/csensitivity[cidx]**4 - 1/ssensitivity[sidx]**4)**(-0.25)
		return sensitivity
	
	def maximum_likelihood_uncertainty_consolidation(self):
		matching_indices = frequency_index_matcher(chunk=chunk)
		axion_fit_uncertainty = self.chunk_axion_fit_uncertainty
		scan_axion_fit_uncertainty = self.scan_axion_fit_uncertainty
		for i, val in enumerate(matching_indices):
			cidx = val[0]
			sidx = val[1]
			if self.op == "+":
				axion_fit_uncertainty[cidx] = (1/axion_fit_uncertainty[cidx]**2 + 1/scan_axion_fit_uncertainty[sidx]**2)**(0.5)
			else:
				axion_fit_uncertainty[cidx] = (1/axion_fit_uncertainty[cidx]**2 - 1/scan_axion_fit_uncertainty[sidx]**2)**(0.5)
		return axion_fit_uncertainty
	
	def model_excess_sqrd_consolidation(self):
		matching_index = frequency_index_matcher(self)
		MES = self.chunk_model_excess_sqrd
		sMES = self.scan_model_excess_sqrd
		
		for i, val in enumerate(matching_index):
			cidx = indices[0]
			sidx = indices[1]
			
			if op=="+":
				MES[cidx] = MES[cidx] + sMES[sidx]
			else:
				MES[cidx] = MES[cidx] - sMES[sidx]
		return MES
	
	def noise_power_consolidation(self):
		SNR = self.chunk_SNR
		WS = self.chunk_optimal_weight_sum
		noise_power = SNR/WS
		return noise_power
	
	def nscan_consolidation(self):
		matching_indices = frequency_index_matcher(chunk=chunk)
		nscan=self.chunk_nscans
		scan_nscan=self.scan_nscans
		for i, val in enumerate(matching_indices):
			cidx = indices[0]
			sidx = indices[1]
			if self.op=="+":
				nscan[cidx] = nscans[cidx] + scan_nscans[sidx]
			else:
				nscans[cidx] = nscans[cidx] - scan_nscans[sidx]
		return nscans
		
	def optimal_weight_sum_consolidation(self):
		matching_indices = frequency_index_matcher(chunk=chunk)
		WS = self.chunk_optimal_weight_sum
		sWS = self.scan_optimal_weight_sum
		for i, val in enumerate(matching_indices):
			cidx = indices[0]
			sidx = indices[1]
			if self.op=="+":
				WS[cidx] = WS[cidx] + sWS[cidx]
			else:
				WS[cidx] = WS[cidx] - sWS[cidx]
		return WS
	
	def power_deviation_consolidation(self):
		matching_indices = frequency_index_matcher(chunk=chunk)
		power_deviation = self.chunk_power_deviation
		spower_deviation = self.scan_power_deviation
		for i, val in enumerate(matching_indices):
			cidx = indices[0]
			sidx = indices[1]
			if self.op=="+":
				power_deviation[cidx] = power_deviation[cidx] + spower_deviation[sidx]
			else:
				power_deviation[cidx] = power_deviation[cidx] - spower_deviation[sidx]
		return power_deviation
	
	def power_senstivity_consolidation(self):
		matching_indices = self.frequency_index_matcher(self)
		sensitivity = self.chunk_sensitivity_power
		ssensitivity = self.scan_sensitivity_power
		for i, val in enumerate(matching_indices):
			cidx = val[0]
			sidx = val[1]
			if self.op=="+":
				sensitivity[cidx]=(1/sensitivity[cidx]**2 + 1/ssensitivity[sidx]**2)**(-0.5)
			else:
				sensitivity[cidx]=(1/sensitivity[cidx]**2 - 1/ssensitivity[sidx]**2)**(-0.5)
		return sensitivity
		
	def sigma_A_consolidation(self):
		matching_indices = frequency_index_matcher(chunk=chunk)
		sigma_A = self.chunk_axion_fit_uncertainty
		ssigma_A = self.scan_axion_fit_uncertainty
		for i, val in enumerate(matching_indices):
			cidx = val[0]
			sidx = val[1]
			if self.op =="+":
				sigma_A[cidx] = (1/(sigma_A[cidx]**2) + 1/(ssigma_A[sidx]**2))**(1/2)
			else:
				sigma_A[cidx] = (1/(sigma_A[cidx]**2) - 1/(ssigma_A[sidx]**2))**(1/2)
		return sigma_A
	
	def weighted_delta_consolidation(self):
		matching_indices = frequency_index_matcher(self)
		weighted_deltas = self.chunk_weighted_deltas
		swds = self.scan_weighted_deltas
		for i, val in enumerate(matching_indices):
			cidx = val[0]
			sidx = val[1]
			if self.op=="+":
				weighted_deltas[cidx] = weighted_deltas[cidx] + swds[sidx]
			else:
				weighted_deltas[cidx] = weighted_deltas[cidx] - swds[sidx]
		return weighted_deltas

		
		
def initialization(chunk):

	if 'scans' not in chunk.attrs:
		chunk.attrs['scans'] = []
		
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
	scan_obj = scan(scan=scan, scid = scan_id, chunk=chunk, op="+")
	
	#Running initiliazation procedure for chunks, i.e. create missing arrays.
	#chunk = initialization(chunk) 									Not currently using.
	
	if op=="+" and scan_id in chunk["scans_in"]:
		phrase = "scan " + str(scan_id) + " already added"
		return phrase
	elif op=="-" and scan_id in chunk["scans_out"]:
		phrase = "scan " + str(scan_id) + " already subtracted"
		return phrase
	
	if type(scan) is str:
		corrupt_scan = True
		add_subtract = 'ommit'
		op = nil
	
	if not corrupt_scan:
		chunk.attrs["nscans"] = scan_obj.nscan_consolidation(chunk,op,scan)
		chunk.attrs["SNR"] = scan_obj.SNR_consolidation(chunk,op,scan)
		chunk.attrs["optimal_weight_sum"] = scan_obj.optimal_weight_sum_consolidation(chunk,op,scan)
		chunk.attrs["noise_power"] = scan_obj.noise_power_consolidation(chunk, op, scan)
		chunk.attrs["model_excess_sqrd"] = scan_obj.model_excess_sqrt_consolidation(chunk, op, scan)
		chunk.attrs["axion_fit_uncertainty"] = scan_obj.sigma_A_consolidation(chunk, op, scan)
		chunk.attrs["axion_fit"] = scan_obj.axion_fit_consolidation(chunk, op, scan)
		chunk.attrs["power_deviation"] = scan_obj.power_devication_consolidation(chunk, op, scan)
		chunk.attrs["axion_fit_significance"] = scan_obj.axion_fit_significance_consolidation(chunk, op, scan)
		chunk.attrs["sensitivity_coupling"] = scan_obj.coupling_sensitivity_consolidation(chunk, op, scan)
		chunk.attrs["sensitivity_power"] = scan_obj.power_senstivity_consolidation(chunk, op, scan)
	
	chunk.attrs["scans"].append(scan_id)
	
	if add_subtract=='add':
		chunk.attrs["scans_id"].append(scan_id)
		for i, chunk_scan_id in enumerate(chunk.attrs["scans_out"]):
			if scan_id == chunk_scan_id:
				chunk.attrs["scans_out"][i]=nil
	elif add_subtract=="subtract":
		chunk.attrs["scans_out"].append(scan_id)
		for i, chunk_scan_id in enumerate(chunk.attrs["scans_in"]):
			if scan_id == chunk_scan_id:
				chunk.attrs["scans_id"][i] = nil
	
	lastcalc = dt.datetime.now()
	lastcalc = lastcalc.strftime('%Y-%m-%d %H:%M:%S')
	chunk["last_change"] = lastcalc
	