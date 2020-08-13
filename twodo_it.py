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

def pseudoit(globs):

	while globs['iteration'] <= globs['num-iters']:
		globs = iterative.mapping(globs);
	# Call the mapping function for each iteration.

	return globs;
#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.
	globs = GV.init();
	
	if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
		sys.exit("# Pesudo-it version " + globs['version'] + " released on " + globs['releasedate']);
	# The version option to simply print the version and exit.

	print("#");
	print("# " + "=" * 125);
	print(PC.welcome());
	if "-h" not in sys.argv:
		print("       Pseudo assembly by iterative mapping.\n");
	# A welcome banner.

	globs = OP.optParse(globs);
	# Getting the input parameters from optParse.

	if globs['norun']:
		sys.exit("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#");

	globs = pseudoit(globs);
	PC.endProg(globs);

#############################################################################
