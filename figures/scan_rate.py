import h5py
import numpy as np
import csv
import matplotlib.pyplot as plt



#Create scan_rate text file

f = h5py.File(b'../data/raw/run1a_data.hdf5', 'r')
dig = f['digitizer_log_run1a']

freq_list = []
n=0
"""
print("Beginning iterations")
for key in dig:
	start = float(dig[key].attrs['start_frequency'])
	stop = float(dig[key].attrs['stop_frequency'])
	res = float(dig[key].attrs['frequency_resolution'])
	attrs = dig[key].attrs
	for freq in np.arange(start, stop, res):
		if not attrs['cut']:
			freq_list.append(freq)
	n+=1
	if n%500==0:
		print(n)
print("Iterations finished. Beginning counting")	
freqs, counts = np.unique(freq_list, return_counts = True)
print("Counting finished. Saving to disk")
np.savetxt('scan_rate.txt', list(zip(freqs, counts)), fmt="%0.12f", delimiter='\t')	

"""
#parse scan_rate file into arrays to plot with
with open('scan_rate.txt', 'r') as f:
	lines = f.readlines()
	for i in range(len(lines)):
		lines[i] = list(map(float, lines[i].split()))
		if lines[:,0]
lines = np.asarray(lines)
plt.title("Scan rate for run 1A")
plt.xlabel("Frequency (MHz)")
plt.ylabel("Number of times frequency scanned")
plt.scatter(lines[:,0], lines[:,1], s=2)
plt.show()
plt.savefig('D:/Users/shaun/Documents/Coding Software/Python/Scripts/New-Analysis-Scheme/oo_analysis/figures/scan_rate.png', dpi = 1000)
