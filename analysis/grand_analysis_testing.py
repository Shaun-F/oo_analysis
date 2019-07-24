




		#begin iteration
		N_iter = len(self.keys)
		scanparam_keys = ['restype', 'notes', 'nbins', 'axion_frequencies_to_inject',
					'pec_vel', 'signal_dataset', 'filter', 'filter_params', 'signal']

		params = dict((key, getattr(self, key)) for key in scanparam_keys)

		#run single analysis routine
		print("\rPerforming single scan analysis       ")
		self.analysis_results = dict((key, analysis(self.dig_dataset[key], scan_number=key, Tsys=self.Tsys[key],**params)) for key in self.keys if not self.dig_dataset[key].attrs['cut'])
		self.analyzed_keys = list(self.analysis_results.keys())
		iterator = range(len(self.analyzed_keys))

		#cycle over scans and make cuts on data

		int_times = dict((key, float(self.dig_dataset[key].attrs['integration_time'])) for key in self.analyzed_keys)
		bandwidths = dict((key, float(self.dig_dataset[key].attrs['frequency_resolution'])) for key in self.analyzed_keys)

		#radiometer dispersion: Tsys/root(B*t) B==bandidth t==integration time
		radiometer_dispersion = dict((key,(1/(int_times[key]*bandwidths[key])**(0.5))) for key in self.analyzed_keys)
		scan_dispersion = dict((key, self.analysis_results[key]['sigma']) for key in self.analyzed_keys)

		toobig_reason = {"True":"Dispersion compared to radiometer dispersion. Too large", "False":""}
		toosmall_reason = {"True": "Dispersion compared to radiometer dispersion. Too small", "False":""}

		disp_toobig = dict((iter, scan_dispersion[iter]>3*radiometer_dispersion[iter]) for iter in self.analyzed_keys)
		disp_toosmall = dict((iter, scan_dispersion[iter]<0.3*radiometer_dispersion[iter]) for iter in self.analyzed_keys)
		cut = dict((key, bool(disp_toobig[key]*disp_toosmall[key])) for key in self.analyzed_keys)
		add_or_subtract= {"True": 'subtract', "False":'add'}

		(add_subtract_scan(add_or_subtract[cut[iter]], self.analysis_results[iter], grand_analysis_class, iter, self.grand_spectra_group, submeta=submeta) for iter in self.analyzed_keys)

		return self.grand_spectra_group, len(np.where(list(cut.values()))[0])













for key in self.keys:
			try:
				#Run single scan analysis
				scan = self.dig_dataset[key]
				if not scan.attrs['cut']:
					scanparam_keys = ['restype', 'notes', 'nbins', 'axion_frequencies_to_inject']
					modparam_keys = ['pec_vel', 'signal_dataset', 'filter', 'filter_params', 'signal']
					
					scanparams = {key: getattr(self, key) for key in scanparam_keys}
					modparams = {key: getattr(self, key) for key in modparam_keys}
					Tsys = {'Tsys': self.Tsys[key]}
					scan_num = {'scan_number':key}
					axion_scan = {'axion_scan': self.dig_dataset[key]}
					
					params = {**scanparams, **modparams, **Tsys, **scan_num, **axion_scan, "submeta":submeta}
					
					#Run single analysis routine
					analysis_start = time.time()
					self.analysis_results[key] =  analysis(scan, **params)
					analysis_stop = time.time()
					
					analysis_timer.append(analysis_stop - analysis_start)
					
					#cycle over scans and make cuts on data
					int_time = float(scan.attrs['integration_time']) #second
					bandwidth = float(scan.attrs["frequency_resolution"])*10**6 #Hz
					
					#Radiometer dispersion: Tsys/root(B*t) B==bandwidth, t==integration time
					radiometer_dispersion = 1/((bandwidth*int_time)**(0.5))
					scan_dispersion = self.analysis_results[key]["sigma"]
					
					cut = False
					cut_reason = ""
					
					if scan_dispersion>3*radiometer_dispersion: #if deltas dispersion is factor of 3 larger than radiometer equation, cut it.
						cut = True
						cut_reason = "Dispersion compared to radiometer dispersion. Too large"
						self.core_analysis.bad_scans.append(key)
					elif scan_dispersion<0.3*radiometer_dispersion:
						cut = True
						cut_reason = "Dispersion compared to radiometer dispersion. Too Small"
						self.core_analysis.bad_scans.append(key)
						print("Background subtraction failed for scan {0}".format(key))
					else:
						pass
						
					#add remaining scans to grand_spectra via coaddition
					coaddition_start = time.time()
					if not scan.attrs['cut']:
						add_subtract_scan('add', self.analysis_results[key], grand_analysis_class, key, self.grand_spectra_group, **{'submeta':submeta})
					
					elif scan.attrs['cut']:
						self.ncut += 1
						add_subtract_scan('subtract', self.analysis_results[key], grand_analysis_class, key, self.grand_spectra_group, **{'submeta':submeta})
					
					if int(counter/N_iter*100)-counter/N_iter*100<10**(-8):
						print("\rPerforming Analysis ( {0}% Complete ) \r".format(int((counter/N_iter)*100)), end='')
					counter+=1
					coaddition_stop = time.time()
					coaddition_timer.append(coaddition_stop - coaddition_start)
					
			except (KeyError, MemoryError, IndexError) as error:
				print("Error at scan {0}. Saving to error log".format(key))
				open('../meta/error_log', 'a+').write(str(time.time())+ "\n\n"+ str(error))
				self.file.close()
				raise
			except KeyboardInterrupt:
				self.file.close()
				print("Interrupted")
				try:
					sys.exit(0)
				except SystemExit:
					os._exit(0)
		iteration_stop = time.time()
		

