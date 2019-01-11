"""
__init__.py: main init file for analysis

Created by: Erik Lentz
Creation Date: 10/26/18
"""

import sys
sys.path.append("..")
import analysis.grand_analysis as grand_analysis


def grand_spectra(object):
	# performs analysis on all scans according to supplied signals
	grand_analysis_class = grand_analysis.analyser(object)
	Grand_Analysis = grand_analysis_class.Grand_Analysis()
	return Grand_Analysis
