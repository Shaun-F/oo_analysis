"""
__init__.py: main init file for analysis

Created by: Erik Lentz
Creation Date: 10/26/18
"""

# undertake analysis calculations
import oo_analysis.control
from oo_analysis.control.core import core_analysis
import time
import argparse
import cProfile, pstats, io, os, sys #Meta analysis

############# Meta analysis
profiler_destination = os.getcwd() + "/oo_analysis/meta/Analytics.txt"

############# Argument parsing
P = argparse.ArgumentParser(description="Main execution file for oo_analysis")
P.add_argument('-t', '--timeit', action='store_true', default=False, help='Time all subprocesses of the analysis \n DEFAULT: False')
P.add_argument('-cgs', '--clear_grand_spectra', action='store_true', default=False, help="Clear the current grand spectra and start from scratch. \n DEFAULT: False")
P.add_argument('--start_scan', action='store', default = '388518', help="Scan number to begin the analysis from If not specified, starting number specified by job.param. \n DEFAULT: 388518")
P.add_argument('--end_scan', action='store', default = '561150', help="Scan number to end the analysis at If not specified, ending number specified by job.param. \n DEFAULT: 561150")
P.add_argument('-p', '--make_plots', action='store_true', default=False, help="Generate plots at end of analysis. \n DEFAULT: False")


args = P.parse_args()
#############



if args.timeit:
	print("################## Running meta analysis ##################")
	profiler = cProfile.Profile()
	profiler.enable()
	
	total_analysis_start = time.time()
	full_analysis = core_analysis(args)

	full_analysis.execute()

	n=1
	while full_analysis.partitioned:
		args.clear_grand_spectra=False
		print("running next partition. (loop number: {0}) (partition [{1}, {2}])".format(n, max(map(int,x.keys))+1, x.end_scan))
		params = {'start_scan': max(map(int,full_analysis.keys))+1, 'end_scan': full_analysis.end_scan}
		full_analysis = core_analysis(**params)
		full_analysis.execute(args)
		n+=1
	full_analysis.garbage_collector()

	total_analysis_stop = time.time()

	print("\nEntire analysis over all {0} partitions took {1:0.3f} seconds".format(n, total_analysis_stop-total_analysis_start))
	
	profiler.disable()
	stream = io.StringIO()
	stats = pstats.Stats(profiler, stream=stream).sort_stats('cumtime')
	stats.print_stats()
	with open(profiler_destination, 'w+') as f:
		f.write(stream.getvalue())
else:
		
	total_analysis_start = time.time()
	full_analysis = core_analysis(args)

	full_analysis.execute()

	n=1
	while full_analysis.partitioned:
		args.clear_grand_spectra=False
		print("running next partition. (loop number: {0}) (partition [{1}, {2}])".format(n, max(map(int,x.keys))+1, x.end_scan))
		params = {'start_scan': max(map(int,full_analysis.keys))+1, 'end_scan': full_analysis.end_scan}
		full_analysis = core_analysis(**params)
		full_analysis.execute(args)
		n+=1
	full_analysis.garbage_collector()

	total_analysis_stop = time.time()

	print("\nEntire analysis over all {0} partitions took {1:0.3f} seconds".format(n, total_analysis_stop-total_analysis_start))
	
