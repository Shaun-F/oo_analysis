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
import cProfile, pstats, io, os #Meta analysis

############# Meta analysis
profiler_destination = os.getcwd() + "/oo_analysis/meta/Analytics.txt"

############# Argument parsing
P = argparse.ArgumentParser(description="Main execution file for oo_analysis")
P.add_argument('-t', '--timeit', action='store_true', default=False, help='Argument specifies whether to time all the subprocesses of the analysis')
P.add_argument('-cgs', '--clear_grand_spectra', action='store_true', default=False, help="Argument specifies whether to delete the grand spectra and start from scratch. Default is False")
P.add_argument('--start_scan', action='store', default = '', help="Argument specifies the starting scan number of the analysis. If not specified, starting number specified by job.param")
P.add_argument('--end_scan', action='store', default = '', help="Argument specifies the ending scan number of the analysis. If not specified, ending number specified by job.param")
P.add_argument('-p', '--make_plots', action='store', default=False, help="Argument specifies whether to generate plots or not after analysis is completed")

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
	stats = pstats.Stats(profiler, stream=stream).sort_stats('time')
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
	
