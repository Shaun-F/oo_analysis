"""
__init__.py: main init file for analysis

Created by: Shaun Fell
"""

# undertake analysis calculations
import oo_analysis.control
from oo_analysis.control.core import core_analysis
import time
import argparse
import cProfile, pstats, io, os, sys #Meta analysis

############# Meta analysis file location
profiler_destination = os.getcwd() + "/oo_analysis/meta/Analytics.txt"

############# Argument parsing
P = argparse.ArgumentParser(description="Main execution file for oo_analysis")
P.add_argument('-t', '--timeit', action='store_true', default=False, help='Time all subprocesses of the analysis \n DEFAULT: False')
P.add_argument('-cgs', '--clear_grand_spectra', action='store_true', default=False, help="Clear the current grand spectra and start from scratch. \n DEFAULT: False")
P.add_argument('--start_scan', action='store', default = '388518', help="Scan number to begin the analysis from If not specified, starting number specified by job.param. \n DEFAULT: 388518")
P.add_argument('--end_scan', action='store', default = '561150', help="Scan number to end the analysis at If not specified, ending number specified by job.param. \n DEFAULT: 561150")
P.add_argument('-p', '--make_plots', action='store_true', default=False, help="Generate plots at end of analysis. \n DEFAULT: False")
P.add_argument('-f', '--filter', action='store', default="RCHPF", help='What background subtraction filter to run the analysis with. \n DEFAULT: RCHPF')

args = P.parse_args()
#############



if args.timeit:
	print("################## Running meta analysis ##################")
	profiler = cProfile.Profile() #Initialize profiler
	profiler.enable() #Enable the profiler
	
	total_analysis_start = time.time() #Get current system time
	full_analysis = core_analysis(args) #Initialize main analysis class

	full_analysis.execute() #Execute the main analysis

	n=1
	while full_analysis.partitioned: #If the data file was partition, run the analysis over the remaining partitions (Paritioning allows the analysis to be run on computers with limited RAM)
		args.clear_grand_spectra=False #Dont clear the grand spectra since were running the analysis on the next partition and will add this partition to the grand spectra
		print("running next partition. (loop number: {0}) (partition [{1}, {2}])".format(n, max(map(int,x.keys))+1, x.end_scan))
		params = {'start_scan': max(map(int,full_analysis.keys))+1, 'end_scan': full_analysis.end_scan} #Next partition
		full_analysis = core_analysis(**params) #Initialize analysis class for next partition
		full_analysis.execute(args) #Execute analysis on next partition
		n+=1
	full_analysis.garbage_collector() #Collect garbage to free up memory

	total_analysis_stop = time.time() #Current system time

	print("\nEntire analysis over all {0} partitions took {1:0.3f} seconds".format(n, total_analysis_stop-total_analysis_start))
	
	profiler.disable() #Disable profiler
	stream = io.StringIO() #Set the input output locations
	stats = pstats.Stats(profiler, stream=stream).sort_stats('cumtime') #Get the profiler statistics
	stats.print_stats() #Print the profile statistics to the input output location
	with open(profiler_destination, 'w+') as f:
		f.write(stream.getvalue()) #Save profiler statistics to file
else:
		
	total_analysis_start = time.time() #Get the current system time
	full_analysis = core_analysis(args) #Initialize the main analysis class

	full_analysis.execute() #Execute the main analysis

	n=1
	#If the data file was partitioned, run the analysis over the remaining partitions (Partitioning allows the analysis to be run on machines with limited RAM)
	while full_analysis.partitioned:
		args.clear_grand_spectra=False #DOnt clear grand spectra since we need to add this next partition into it
		print("running next partition. (loop number: {0}) (partition [{1}, {2}])".format(n, max(map(int,x.keys))+1, x.end_scan))
		params = {'start_scan': max(map(int,full_analysis.keys))+1, 'end_scan': full_analysis.end_scan} #Next partition
		full_analysis = core_analysis(**params) #Initialize main analysis class for next partition
		full_analysis.execute(args) #Execute the main analysis on this partition
		n+=1
	full_analysis.garbage_collector() #Collect garbage to free up system memory

	total_analysis_stop = time.time() #Current system time

	print("\nEntire analysis over all {0} partitions took {1:0.3f} seconds".format(n, total_analysis_stop-total_analysis_start))
	
