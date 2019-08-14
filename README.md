# Object-Oriented Analysis

Object-Oriented Analysis (OO-Analysis) is a python package for analyzing ADMX data offline.

## Installation

Clone the repository from Github (https://github.com/Shaun-F/oo_analysis). The data files are too large to be stored on Github, so contact me directly and we can figure out a way to send them to you.
Please see the oo_analysis/meta/requirements.txt file for a complete list of currently implemented packages and required programs.

## Usage

OO_Analysis can be run directly from the command line equipped a python 3 interpreter. 
```bash
python -m oo_analysis
```

OO_Analysis can also be run from within the python 3 interpreter just like any other python package. Arguments have to passed to the class during initialization.

```python
#First import the package
import oo_analysis

#Then initialize the main execution class
cl = oo_analysis.core(make_plots=True)

#Then execute the analysis
cl.execute()
```

There are a number of arguments that can be passed to the package that dictates the execution. 
To see the full list of available arguments, run
```bash
python -m oo_analysis -h
```

## Contributing

Pull requests are welcome. 

