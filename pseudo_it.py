#!/usr/bin/env python3
#############################################################################
# Iterative read-mapping and variant calling for pseudo-reference assembly
# This is the main interface.
#
# Gregg Thomas
# Fall 2019
#############################################################################

import sys
import os
import multiprocessing as mp
import shutil
import pseudo_it_lib.picore as PC
import pseudo_it_lib.opt_parse as OP
import pseudo_it_lib.global_vars as GV
import pseudo_it_lib.pimap as pimap
import pseudo_it_lib.piref as piref
import pseudo_it_lib.iterative as iterative

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
		print("# Pesudo-it version " + globs['version'] + " released on " + globs['releasedate'])
		sys.exit(0);
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
		print("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#")
		sys.exit(0);

	globs = pseudoit(globs);
	PC.endProg(globs);

#############################################################################
