"""
param_parser.py: module for parsing parameter file

Created by: Erik Lentz
Creation Date: 10/26/18
"""

def parser(filename):
    # function takes filename and parses it into dict for analysis parameter setting

    # open, read file line by line, adding params to dictinary
    file = open(filename,'r')
    params = {}
    for line in file:
        if '#' not in line and line!="\n":
            split = line.split('=')
            params[split[0]] = eval(split[1][:-1]) #-1 accounts for return at end of line. eval() converts string to python data structure

    # close file
    file.close()
    return params
