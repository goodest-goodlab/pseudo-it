# Parsing and printing the options and meta-info for Pseudo-it.
# Much of the error checking is done here as well.
#############################################################################
import sys, os, math, lib.picore as PC
#############################################################################

def optParse(globs):
# This function handles the command line options and prepares the output directory and files.
	try:
		import argparse;
	except:
		PC.errorOut("\n*** ERROR: Your installation of Python is missing the argparse module. Please try a different version of Python (3+) or install the module.\n")
	# First check if the argparse module is installed. If not, the input options cannot be parsed.
	try:
		import psutil
		globs['psutil'] = True;
	except:
		globs['psutil'] = False;
	# Check if psutil is installed for memory usage stats.

	parser = argparse.ArgumentParser(description="Pseudo-it: assembly by iterative mapping.");

	parser.add_argument("-ref", dest="ref", help="The FASTA assembly to use for the initial mapping.", default=False);
	parser.add_argument("-se", dest="se", help="A FASTQ file containing single-end reads.", default=False);
	parser.add_argument("-pe1", dest="pe1", help="A FASTQ file containing pair 1 of paired-end reads.", default=False);
	parser.add_argument("-pe2", dest="pe2", help="A FASTQ file containing pair 2 of paired-end reads.", default=False);
	parser.add_argument("-pem", dest="pem", help="A FASTQ file containing merged paired-end reads.", default=False);
	# Inputs
	parser.add_argument("-tmp", dest="tmp_dir", help="Some programs write files to a temporary directory. If your default tmp dir is size limited, specify a new one here. Default: Your system's tmp dir.", default=False);
	parser.add_argument("-o", dest="out_dest", help="Desired output directory. Default: pseudoit-[date]-[time]", default=False);
	# Output
	parser.add_argument("-resume", dest="resume", help="The path to a previous Pseudo-it directory to resume a run. Scans for presence of files and resumes when it can't find an expected file.", default=False);
	parser.add_argument("-bwa", dest="bwa_path", help="The path to the BWA mapping progam. Default: bwa", default=False);
	parser.add_argument("-picard", dest="picard_path", help="The exact command used to run picard. For a jar file: java -jar <full path to jar file>. For an alias or conda install: picard. Include heap size in command, i.e. -Xmx6g. Default: picard", default=False);
	parser.add_argument("-samtools", dest="samtools_path", help="The path to the samtools progam. Default: samtools", default=False);
	parser.add_argument("-gatk", dest="gatk_path", help="The path to the GATK progam. Default: gatk", default=False);
	parser.add_argument("-bedtools", dest="bedtools_path", help="The path to the bedtools progam. Default: bedtools", default=False);
	parser.add_argument("-bcftools", dest="bcftools_path", help="The path to the bcftools progam. Default: bcftools", default=False);
	# Dependency paths
	parser.add_argument("-i", dest="num_iters", help="The number of iterations Pseudo-it will run. Default: 4.", default=False);
	parser.add_argument("-bwa-t", dest="bwa_threads", help="The number of threads for BWA mem to use for each library. If you specify -bwa-t 3 and have 3 libraries (and have at least -p 3), this means a total of 9 processes will be used during mapping. If left unspecified and -p is specified this will be determined automatically by dividing -p by the number of libraries you provide. Otherwise, default: 1.", default=False);
	parser.add_argument("-gatk-t", dest="gatk_threads", help="The number of threads for GATK's Haplotype caller to use. If you specify -p 4 and -gatk-t 4, this means that a total of 16 processes will be used. GATK default: 4.", default=False);
	parser.add_argument("-heap", dest="java_heap", help="The heap size to allot for called java programs (picard). Enter an integer, assumed unit is gigabytes (g). Default: System default.", default=False);
	parser.add_argument("-f", dest="filter", help="The expression to filter variants. Must conform to VCF INFO field standards. Default: \"MQ < 30.0 || DP < 5 || DP > 60\"", default=False);
	parser.add_argument("-p", dest="processes", help="The number of processes Pseudo-it should use. Default: 1.", default=False);
	# User params
	parser.add_argument("--noindels", dest="noindels", help="Set this to not incorporate indels into the final assembly.", action="store_true", default=False);
	parser.add_argument("--overwrite", dest="ow_flag", help="Set this to overwrite existing files.", action="store_true", default=False);
	parser.add_argument("--quiet", dest="quiet_flag", help="Set this flag to prevent Referee from reporting detailed information about each step.", action="store_true", default=False);
	parser.add_argument("--version", dest="version_flag", help="Simply print the version and exit. Can also be called as '-version', '-v', or '--v'", action="store_true", default=False);
	# User options
	parser.add_argument("--norun", dest="norun", help=argparse.SUPPRESS, action="store_true", default=False);
	parser.add_argument("--dryrun", dest="dryrun", help=argparse.SUPPRESS, action="store_true", default=False);
	parser.add_argument("--debug", dest="debug_opt", help=argparse.SUPPRESS, action="store_true", default=False);
	parser.add_argument("--nolog", dest="nolog_opt", help=argparse.SUPPRESS, action="store_true", default=False);
	# Performance tests
	args = parser.parse_args();
	# The input options and help messages

	if args.norun:
		globs['norun'] = True;
		globs['log-v'] = -1;
	if args.dryrun:
		globs['dryrun'] = True;
	globs['overwrite'] = args.ow_flag;

	if not any([args.se, args.pe1, args.pe2, args.pem]):
		PC.errorOut("OP1", "At least one FASTQ library file must be given with -se, -pe1 and -pe2, or -pem", globs);

	if args.pe1 and not args.pe2 or not args.pe1 and args.pe2:
		PC.errorOut("OP2", "With a paired end library, both -pe1 and -pe2 must be specified.", globs);

	if not args.ref:
		PC.errorOut("OP3", "A reference FASTA file must be provided with -ref.", globs);

	globs['se'], globs['pe1'], globs['pe2'], globs['pem'], globs['ref'] = args.se, args.pe1, args.pe2, args.pem, args.ref;
	globs['num-libs'] = len( [ l for l in [globs['se'], globs['pe1'], globs['pem']] if l ] );

	for l in ['se', 'pe', 'pem']:
		if l != 'pe':
			if globs[l]:
				globs['libs'][l] = globs[l];
		elif globs[l+"1"]:
			globs['libs'][l] = globs[l+"1"] + " " + globs[l+"2"];
	# Restructure the FASTQ libraries into a single dictionary.

	PC.fileCheck(globs);
	# Input file checking.

	globs = PC.execCheck(globs, args);
	# Check the dependency paths.

	if args.resume:
		if globs['overwrite']:
			PC.errorOut("OP4", "--overwrite cannot be specified with --resume.", globs);
		if not os.path.isdir(args.resume):
			PC.errorOut("OP5", "Directory specified by --resume does not exist!", globs);
		globs['resume'] = True;
		globs['outdir'] = args.resume;
	# Check for resume flag.

	else:
		if not args.out_dest:
			globs['outdir'] = "pseudoit-out-" + globs['startdatetime'];
		else:
			globs['outdir'] = args.out_dest;

		if not globs['overwrite'] and os.path.exists(globs['outdir']):
			PC.errorOut("OP6", "Output directory already exists: " + globs['outdir'] + ". Specify new directory name, set --overwrite to overwrite, or set -resume to resume.", globs);

		if not os.path.isdir(globs['outdir']) and not globs['norun']:
			os.makedirs(globs['outdir']);

	globs['logdir'] = os.path.join(globs['outdir'], "logs");
	if not os.path.isdir(globs['logdir']) and not globs['norun']:
		os.makedirs(globs['logdir']);

	if args.tmp_dir:
		globs['tmpdir'] = args.tmp_dir;
		if not os.path.isdir(globs['tmpdir']):
			os.system("mkdir " + globs['tmpdir']);

	globs['logfilename'] = os.path.join(globs['outdir'], os.path.basename(os.path.normpath(globs['outdir'])) + ".log");
	if globs['dryrun']:
		globs['logfilename'] = globs['logfilename'].replace(".log", "-dryrun.log");
	globs['scaffs'] = os.path.join(globs['outdir'], "scaffold-list.txt");
	globs['endprog'] = True;
	# Output prep.

	if args.noindels:
		globs['indels'] = False;
	# Indel output option.

	if args.debug_opt:
		globs['debug'] = True;
	if args.nolog_opt:
		globs['log-v'] = -1;
	# Hidden test options

	if args.num_iters and ("-" in args.num_iters or not args.num_iters.isdigit()):
		PC.errorOut("OP7", "-i must be an integer value greater than 1.", globs);
	elif args.num_iters:
		globs['num-iters'] = int(args.num_iters);
	# Checking the number of iterations option.

	# if args.java_heap and not isdigit(args.java_heap) or int(args.java_heap) < 1:
	# 	PC.errorOut("OP8", "-heap must be a positive integer.", globs);
	# elif args.java_heap:
	# 	globs['heap'] = int(args.java_heap);
	# Getting the java heap option.

	if args.filter:
		globs['filter'] = args.filter;
	# Check the filter option.

	if args.processes and ("-" in args.processes or not args.processes.isdigit()):
		PC.errorOut("OP8", "-p must be an integer value greater than 1.", globs);
	elif args.processes:
		globs['num-procs'] = int(args.processes);
	# Checking the number of processors option.

	if args.bwa_threads:
		if "-" in args.bwa_threads or not args.bwa_threads.isdigit():
			PC.errorOut("OP9", "-bwa-t must be an integer value greater than 1.", globs);
		else:
			globs['bwa-t'] = int(args.bwa_threads);
	elif args.processes != 1:
		globs['bwa-t'] = math.floor(globs['num-procs'] / globs['num-libs']);
	# Getting the number of BWA mem threads.

	if args.gatk_threads:
		if "-" in args.gatk_threads or not args.gatk_threads.isdigit():
			PC.errorOut("OP9", "-gatk-t must be an integer value greater than 1.", globs);
		else:
			globs['gatk-t'] = int(args.gatk_threads);

	if args.quiet_flag:
		globs['quiet'] = True;
	# Check the quiet option

	startProg(globs);
	# After all the essential options have been set, call the welcome function.

	if args.quiet_flag:
		globs['log-v'] = 3;
	# If the quiet option was set before, set the verbosity for printWrite here so nothing is
	# printed to the screen, but still output to the logfile.

	step_start_time = "";
	if globs['psutil']:
		globs['pids'] = [psutil.Process(os.getpid())];	
	globs['stats'] = True;
	if not globs['norun']:
		step_start_time = PC.report_stats(globs, stat_start=True);
	# Initializing the stats options if --quiet is not set.

	return globs, step_start_time;

#############################################################################

def startProg(globs):
# A nice way to start the program.
	print("#");
	PC.printWrite(globs['logfilename'], globs['log-v'], "# Welcome to Pseudo-it -- Pseudo genome assembly via iterative mapping.");
	PC.printWrite(globs['logfilename'], globs['log-v'], "# Version " + globs['version'] + " released on " + globs['releasedate']);
	PC.printWrite(globs['logfilename'], globs['log-v'], "# Pseudo-it was developed by Brice Sarver, Gregg Thomas, and Jeffrey Good");
	PC.printWrite(globs['logfilename'], globs['log-v'], "# Citation:      " + globs['doi']);
	PC.printWrite(globs['logfilename'], globs['log-v'], "# Website:       " + globs['http']);
	PC.printWrite(globs['logfilename'], globs['log-v'], "# Report issues: " + globs['github']);
	PC.printWrite(globs['logfilename'], globs['log-v'], "#");
	PC.printWrite(globs['logfilename'], globs['log-v'], "# The date and time at the start is: " + PC.getDateTime());
	PC.printWrite(globs['logfilename'], globs['log-v'], "# Using Python version:              " + globs['pyver'] + "\n#");
	PC.printWrite(globs['logfilename'], globs['log-v'], "# The program was called as:         " + " ".join(sys.argv) + "\n#");

	pad = 35;
	opt_pad = 50;
	PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
	PC.printWrite(globs['logfilename'], globs['log-v'], "# INPUT/OUTPUT INFO:");
	
	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Reference FASTA file:", pad) + globs['ref']);
	
	if globs['se']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Single-end reads:", pad) + globs['se']);
	if globs['pem']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Merged paired-reads:", pad) + globs['pem']);	
		if globs['pe1']:
			PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Unmerged paired-reads 1:", pad) + globs['pe1']);	
			PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Unmerged paired-reads 2:", pad) + globs['pe2']);
	elif globs['pe1']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Paired-reads 1:", pad) + globs['pe1']);	
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Paired-reads 2:", pad) + globs['pe2']);

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Temporary file directory:", pad) + globs['tmpdir']);
	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Output directory:", pad) + globs['outdir']);
	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Log file:", pad) + os.path.basename(globs['logfilename']));

	PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
	PC.printWrite(globs['logfilename'], globs['log-v'], "# DEPENDENCY PATHS:");	
	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Program", pad) + "Specified Path");
	# PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
	
	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -bwa", pad) + globs['bwa-path']);
	# Reporting BWA path.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -picard", pad) + globs['picard-path']);
	# Reporting Picard path.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -samtools", pad) + globs['samtools-path']);
	# Reporting samtools path.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -gatk", pad) + globs['gatk-path']);
	# Reporting GATK path.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -bedtools", pad) + globs['bedtools-path']);
	# Reporting bedtools path.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -bcftools", pad) + globs['bcftools-path']);
	# Reporting bcftools path.

	PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
	PC.printWrite(globs['logfilename'], globs['log-v'], "# OPTIONS INFO:");	
	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Option", pad) + PC.spacedOut("Current setting", opt_pad) + "Current action");
	# PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);

	if globs['resume']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -resume", pad) +
					PC.spacedOut(str(globs['outdir']), opt_pad) + 
					"Pseudo-it will attempt to resume the run from this directory.");
	# Reporting the resume option.

	if globs['indels']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --noindels", pad) + 
					PC.spacedOut("False", opt_pad) + 
					"Final assembly will incorporate indels.");
	else:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --noindels", pad) + 
					PC.spacedOut("True", opt_pad) + 
					"Final assembly will NOT incorporate indels.");		
	# Reporting --noindels option.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -i", pad) + 
				PC.spacedOut(str(globs['num-iters']), opt_pad) + 
				"Pseudo-it will perform this many iterations of mapping.");
	# Reporting the number of iterations.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -f", pad) + 
				PC.spacedOut(str(globs['filter']), opt_pad) + 
				"Variants will be filtered based on these criteria during the Variant Filtration step.");
	# Reporting the filter option.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -p", pad) + 
				PC.spacedOut(str(globs['num-procs']), opt_pad) + 
				"Pseudo-it will use this many processes to run.");
	# Reporting the processes option.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -bwa-t", pad) + 
				PC.spacedOut(str(globs['bwa-t']), opt_pad) + 
				"BWA mem will use this many threads for mapping.");
	# Reporting the BWA mem threads option.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -gatk-t", pad) + 
				PC.spacedOut(str(globs['gatk-t']), opt_pad) + 
				"HaplotypeCaller's --native-pair-hmm-threads option will use this many threads.");
	# Reporting the GATK HaplotypeCaller --native-pair-hmm-threads option.

	# PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -heap", pad) + 
	# 			PC.spacedOut(str(globs['heap']), pad) + 
	# 			"The heap size in gigabytes for called java programs (picard).");
	# Reporting the GATK HaplotypeCaller --native-pair-hmm-threads option.

	if not globs['quiet']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --quiet", pad) + 
					PC.spacedOut("False", opt_pad) + 
					"Runtime, memory, and command info will be printed to the screen while Pseudo-it is running.");
	else:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --quiet", pad) + 
					PC.spacedOut("True", opt_pad) + 
					"No further information will be printed to the screen while Pseudo-it is running.");
		PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
		PC.printWrite(globs['logfilename'], globs['log-v'], "# Running...");
	# Reporting the quiet option.

	if globs['debug']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --debug", pad) + 
					PC.spacedOut("True", opt_pad) + 
					"Printing out a bit of debug info.");
	# Reporting the debug option.

	if globs['norun']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --norun", pad) + 
					PC.spacedOut("True", opt_pad) + 
					"ONLY PRINTING RUNTIME INFO.");
		PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
	# Reporting the norun option.

#############################################################################