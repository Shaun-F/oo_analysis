"""
scan_analysis.py: manages the analysis of all scan_analysis

Created by: Erik Lentz
"""
import oo_analysis.analysis.MR_scan_analysis
from oo_analysis.analysis.MR_scan_analysis import MR_scan_analysis

# handle the background-subtracted data\

def analysis(scan, **params):
	if params['restype'] == "HR":
		pass
		#Run high resolution analysis routine
		
	if params['restype'] == "MR":
		return MR_scan_analysis(scan, **params)
		
	if params['restype'] == "LR":
		pass
		#Run low resolution analysis routine
	
