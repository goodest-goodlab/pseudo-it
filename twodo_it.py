#!/usr/bin/python
#############################################################################
# Iterative read-mapping and variant calling for pseudo-reference assembly
# This is the main interface.
#
# Gregg Thomas
# Fall 2019
#############################################################################

import sys, os, multiprocessing as mp, shutil, lib.picore as PC, \
	lib.opt_parse as OP, lib.global_vars as GV, lib.pimap as pimap, \
	lib.piref as piref, lib.iterative as iterative

#############################################################################

def pseudoit(globs, step_start_time):

	if globs['psutil']:
		import psutil
	step_start_time = PC.report_stats(globs, "Starting", step_start=step_start_time);
	# Initialize the stats output if --stats is set.

	while globs['iteration'] <= globs['num-iters']:
		globs, step_start_time = iterative.mapping(globs, step_start_time);
	# Call the mapping function for each iteration.

	if globs['stats']:
		step_start_time = PC.report_stats(globs, "End program", step_start=step_start_time, stat_end=True);
	# A step update for --stats.

	return;
#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.
	globs = GV.init();
	
	if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
		sys.exit("# Pesudo-it version " + globs['version'] + " released on " + globs['releasedate']);
	# The version option to simply print the version and exit.

	print("#");
	print("# " + "=" * 100);
	print(PC.welcome());
	print("       Pseudo assembly by iterative mapping.\n")
	# A welcome banner.

	globs, step_start_time = OP.optParse(globs);
	# Getting the input parameters from optParse.

	if globs['norun']:
		sys.exit("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#");

	pseudoit(globs, step_start_time);
	PC.endProg(globs);

#############################################################################
