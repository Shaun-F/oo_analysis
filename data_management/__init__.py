"""
__init__.py: main init file for data management

Created by: Erik Lentz
Creation Date: 10/26/18
"""
import h5py

def write_out(dataset,path):
    data = h5py.File(path, 'w')
    # use write functioins to put grand spectra to file
    return None
