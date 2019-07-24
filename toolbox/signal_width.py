import numpy as np
import numba

def calc_signal_width(signal):
	max = np.max(signal)
	min = np.min(signal)
	half_max = (max-min)/2 + min
	subtracted = signal - half_max
	intersects = np.where(np.diff(np.sign(subtracted)))[0]
	if len(intersects)==1:
		return intersects[0]
	else:
		return intersects[1] - intersects[0]