import h5py
import numpy as np
import csv
import matplotlib.pyplot as plt



#Create scan_rate text file

f = h5py.File(b'../data/raw/run1a_data.hdf5', 'r')
dig = f['digitizer_log_run1a']

freq_list = []
analysis_stop = 680*10**6
analysis_start = 645*10**6
analysis_res = 95.4
domain = np.arange(analysis_start, analysis_stop, analysis_res)
counter = [0]*len(domain)
N = int(((analysis_stop-analysis_start)/analysis_res))
print(N)
print("Beginning iterations")
for key in dig:
	start = float(dig[key].attrs['start_frequency'])*10**6
	stop = float(dig[key].attrs['stop_frequency'])*10**6
	res = float(dig[key].attrs['frequency_resolution'])*10**6
	attrs = dig[key].attrs
	try:
		if not attrs['cut']:
			for freq in np.arange(start, stop, res):
				inx = int((freq-analysis_start)/analysis_res)
				if inx>=N or inx<0:
					print(inx)
					pass
				else:
					counter[inx]+=1
	except KeyError:
		pass
print("Iterations finished. Beginning counting")	


#parse scan_rate file into arrays to plot with

plt.title("Scan rate for run 1A")
plt.xlabel("Frequency (MHz)")
plt.ylabel("Number of times frequency scanned")
plt.scatter(domain, counter, s=2)
plt.show()
plt.savefig('D:/Users/shaun/Documents/Coding Software/Python/Scripts/New-Analysis-Scheme/oo_analysis/figures/scan_rate.pdf', dpi = 1000)
