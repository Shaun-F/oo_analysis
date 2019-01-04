import numpy as np

def axion_power(ch, scan_number, form_factor, axion_coupling, halo_density):
	if ch!=2:
		ch=1
	if not form_factor:
		if ch==1:
			form_factor=0.4
		else:
			form_factor=0.08
	
	if not axion_coupling:
		axion_coupling=0.36 #KSVZ is 0.97 DFSZ is 0.36
	if not halo_density:
		halo_density = 0.45 # GeV c^-2 cm^-3
	
	#get axion dataset
	file = h5py.File('run1a_data', 'r')
	scan = file["axion_log_run1a"][str(scan_number)]
	if scan=="error":
		return "error: Fault in axion dataset fetching"
	
	#extract relevent experimental parameters
	q = scan.attrs["Q"]
	B = scan.attrs["b_field"]
	vol = 220 
	res_freq = scan.attrs["mode_frequency"]
	
	total_power = (2.09*10**(-22))*(vol/220)*((B/7.6)**2)*form_factor*((axion_coupling/0.36)**2)*(halo_density/0.45)*(res_freq/750)*(q/70000)*f
	
	
	