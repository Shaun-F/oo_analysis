"""
__init__.py: main init file for analysis

Created by: Erik Lentz
Creation Date: 10/26/18
"""
from .grand_analysis import analyzer


def grand_spectra(object):
	# performs analysis on all scans according to supplied signals
	grand_analysis_class = analyzer(object)
	Grand_Analysis = grand_analysis_class.Grand_Analysis()
	return Grand_Analysis
