"""
__init__.py: main init file for analysis

Created by: Erik Lentz
Creation Date: 10/26/18
"""
import grand_analysis

def grand_spectra(dig_dataset,axion_dataset,signal_dataset):
    # performs analysis on all scans according to supplied signals
    return grand_analysis.grand_spectra(dig_dataset,axion_dataset,signal_dataset)
