# Parsing and printing the options and meta-info for Pseudo-it.
# Much of the error checking is done here as well.
#############################################################################
import sys, os, argparse, math, lib.picore as PC
#############################################################################

def optParse(globs):
# This function handles the command line options and prepares the output directory and files.
	parser = argparse.ArgumentParser(description="Pseudo-it: Genomes assembly by iterative mapping.");

	parser.add_argument("-ref", dest="ref", help="The FASTA assembly to use for the initial mapping.", default=False);
	parser.add_argument("-se", dest="se", help="A FASTQ file containing single-end reads. At least one of -se, -pe1 and pe2, or pem must be provided.", default=False);
	parser.add_argument("-pe1", dest="pe1", help="A FASTQ file containing pair 1 of paired-end reads. At least one of -se, -pe1 and pe2, or pem must be provided.", default=False);
	parser.add_argument("-pe2", dest="pe2", help="A FASTQ file containing pair 2 of paired-end reads. At least one of -se, -pe1 and pe2, or pem must be provided.", default=False);
	parser.add_argument("-pem", dest="pem", help="A FASTQ file containing merged paired-end reads. At least one of -se, -pe1 and pe2, or pem must be provided.", default=False);
	parser.add_argument("-bam", dest="bam", help="OPTIONAL: A BAM file with the provided reads mapped to the reference to be used for the first iteration. This BAM file must be pre-indexed with samtools index.", default=False);
	# Inputs
	parser.add_argument("-tmp", dest="tmp_dir", help="Some programs write files to a temporary directory. If your default tmp dir is size limited, specify a new one here, or just specifiy 'tmp-pi-out' to have a folder called 'tmp' created and used within the main output folder. Default: Your system's tmp dir.", default=False);
	parser.add_argument("-o", dest="out_dest", help="Desired output directory. This will be created for you if it doesn't exist. Default: pseudoit-[date]-[time]", default=False);
	# Output
	parser.add_argument("-resume", dest="resume", help="The path to a previous Pseudo-it directory to resume a run. Scans for presence of files and resumes when it can't find an expected file.", default=False);
	parser.add_argument("-mapper", dest="mapper", help="The name of the mapping progam. One of: bwa, hisat2. Default: bwa", default=False);
	parser.add_argument("-mapperpath", dest="mapper_path", help="The path to the mapping program. By default, just the name of the program.", default=False);
	#parser.add_argument("-bwa", dest="bwa_path", help="The path to the BWA mapping progam. Default: bwa", default=False);
	parser.add_argument("-picard", dest="picard_path", help="The exact command used to run picard. For a jar file: java -jar <full path to jar file>. For an alias or conda install: picard. Include heap size in command, i.e. -Xmx6g. Default: picard", default=False);
	parser.add_argument("-samtools", dest="samtools_path", help="The path to the samtools progam. Default: samtools", default=False);
	parser.add_argument("-gatk", dest="gatk_path", help="The path to the GATK progam. Default: gatk", default=False);
	parser.add_argument("-bedtools", dest="bedtools_path", help="The path to the bedtools progam. Default: bedtools", default=False);
	parser.add_argument("-bcftools", dest="bcftools_path", help="The path to the bcftools progam. Default: bcftools", default=False);
	# Dependency paths
	parser.add_argument("-i", dest="num_iters", help="The number of iterations Pseudo-it will run. Default: 4.", type=int, default=4);
	parser.add_argument("-map-t", dest="map_threads", help="The number of threads for the mapper to use for each read type. If you specify -map-t 3 and have 3 read types (and have at least -p 3), this means a total of 9 processes will be used during mapping. If left unspecified and -p is specified this will be determined automatically by dividing -p by the number of libraries you provide. Otherwise, default: 1.", default=False);
	parser.add_argument("-gatk-t", dest="gatk_threads", help="The number of threads for GATK's Haplotype caller to use. If you specify -p 4 and -gatk-t 4, this means that a total of 16 processes will be used. GATK default: 4.", default=False);
	parser.add_argument("-vcf", dest="vcf", help="A VCF with variants to IGNORE in the called data.", default=False);
	parser.add_argument("-f", dest="filter", help="The expression to filter variants. Must conform to VCF INFO field standards. Default read depth filters are optimized for a 30-40X sequencing run -- adjust for your assembly. Default: \"MQ < 30.0 || DP < 5 || DP > 60\"", default=False);
	parser.add_argument("-p", dest="processes", help="The MAX number of processes Pseudo-it can use. If -p is set to 12 and -gatk-t is set to 4, then Pseudo-it will spawn 3 GATK processes in parallel. Default: 1.", type=int, default=1);
	# User params
	parser.add_argument("--nomkdups", dest="no_mkdups", help="Do not run Picard's MarkDuplicates on mapped reads.", action="store_true", default=False);
	parser.add_argument("--hardmask", dest="hard_mask", help="By default, low quality positions in the final consensus sequence will be soft-masked as lower-case nucleotides. Set this to hard-mask them as Ns.", action="store_true", default=False);
	parser.add_argument("--maponly", dest="map_only", help="Only do one iteration and stop after read mapping.", action="store_true", default=False);
	parser.add_argument("--noindels", dest="noindels", help="Set this to not incorporate indels into the final assembly.", action="store_true", default=False);
	parser.add_argument("--diploid", dest="diploid", help="Set this use IUPAC ambiguity codes in the final FASTA file.", action="store_true", default=False);
	parser.add_argument("--keepall", dest="keep_all", help="By default, pseudo-it keeps only the final files for each step of each iteration (BAM, VCF, FASTA and their respective indices). Set this option to keep all intermediate files. While this is the best way to ensure your runs can be resumed with different settings this will result many large files being saved (total of ~1TB for a 30X genome and 4 iterations).", action="store_true", default=False);
	parser.add_argument("--keeponlyfinal", dest="keep_only_final", help="By default, pseudo-it keeps only the final files for each step of each iteration (BAM, VCF, FASTA and their respective indices). Set this option to keep these files ONLY for the final iteration. While this minimizes storage space required, you will be unable to resume this run.", action="store_true", default=False);
	parser.add_argument("--overwrite", dest="ow_flag", help="Set this to overwrite existing files.", action="store_true", default=False);
	parser.add_argument("--quiet", dest="quiet_flag", help="Set this flag to prevent psuedo-it from reporting detailed information about each step.", action="store_true", default=False);
	parser.add_argument("--version", dest="version_flag", help="Simply print the version and exit. Can also be called as '-version', '-v', or '--v'", action="store_true", default=False);
	# User options
	parser.add_argument("--depcheck", dest="depcheck", help="Run this to check that all dependencies are installed at the provided path. No other options necessary.", action="store_true", default=False);
	parser.add_argument("--dryrun", dest="dryrun", help="With all options provided, set this to run through the whole pseudo-it pipeline without executing external commands.", action="store_true", default=False);
	# Run options
	parser.add_argument("--norun", dest="norun", help=argparse.SUPPRESS, action="store_true", default=False);
	# parser.add_argument("--debug", dest="debug_opt", help=argparse.SUPPRESS, action="store_true", default=False);
	parser.add_argument("--nolog", dest="nolog_opt", help=argparse.SUPPRESS, action="store_true", default=False);
	# Performance tests
	args = parser.parse_args();
	# The input options and help messages

	globs, deps_passed = PC.execCheck(globs, args);
	if args.depcheck:
		if deps_passed:
			print("\n# All dependencies PASSED.\n")
			sys.exit(0);
		else:
			print("\n# Some dependencies NOT FOUND. Please check your installations and provided paths.\n");
			sys.exit(1);
	# Check the dependency paths.		

	if args.norun:
		globs['norun'] = True;
		globs['log-v'] = -1;
	if args.dryrun:
		globs['dryrun'] = True;
	globs['overwrite'] = args.ow_flag;
	# Check run mode options.

	if args.num_iters > 1 or not args.bam:
		if not any([args.se, args.pe1, args.pe2, args.pem]):
			PC.errorOut("OP1", "At least one FASTQ library file must be given with -se, -pe1 and -pe2, or -pem", globs);

	if args.pe1 and not args.pe2 or not args.pe1 and args.pe2:
		PC.errorOut("OP2", "With a paired end library, both -pe1 and -pe2 must be specified.", globs);

	if not args.ref:
		PC.errorOut("OP3", "A reference FASTA file must be provided with -ref.", globs);
	# Check that all input file options are provided.

	if args.bam:
		globs['bam'] = args.bam;
		globs['bam-index'] = args.bam + ".bai";
	# Options if a pre-mapped BAM file is supplied

	globs['se'], globs['pe1'], globs['pe2'], globs['pem'], globs['ref'] = args.se, args.pe1, args.pe2, args.pem, args.ref;
	globs['num-libs'] = len( [ l for l in [globs['se'], globs['pe1'], globs['pem']] if l ] );
	# Parse the fastq files.

	globs = PC.fileCheck(globs);
	# Input file checking.

	if args.no_mkdups:
		globs['mkdups'] = False;
	if args.hard_mask:
		globs['softmask'] = False;
	if args.vcf:
		if not os.path.isfile(args.vcf):
			PC.errorOut("OP4", "Cannot find file specified with -vcf: " + args.vcf, globs);
		globs['in-vcf'] = args.vcf;
	# Options for marking duplicates, hard vs. soft masking, and a vcf file to filter SNPs with.

	for l in ['se', 'pe', 'pem']:
		if l != 'pe':
			if globs[l]:
				globs['libs'][l] = globs[l];
		elif globs[l+"1"]:
			globs['libs'][l] = globs[l+"1"] + " " + globs[l+"2"];
	# Restructure the FASTQ libraries into a single dictionary.

	if args.resume:
		if globs['overwrite']:
			PC.errorOut("OP5", "--overwrite cannot be specified with --resume.", globs);
		if not os.path.isdir(args.resume):
			PC.errorOut("OP6", "Directory specified by --resume does not exist!", globs);
		globs['resume'] = True;
		globs['outdir'] = args.resume;
	# Check for resume flag.

	else:
		if not args.out_dest:
			globs['outdir'] = "pseudoit-out-" + globs['startdatetime'];
		else:
			globs['outdir'] = args.out_dest;

		if not globs['overwrite'] and os.path.exists(globs['outdir']):
			PC.errorOut("OP7", "Output directory already exists: " + globs['outdir'] + ". Specify new directory name, set --overwrite to overwrite, or set -resume to resume.", globs);

		if not os.path.isdir(globs['outdir']) and not globs['norun']:
			os.makedirs(globs['outdir']);
	# Main output dir

	if args.tmp_dir:
		if args.tmp_dir == "tmp-pi-out":
			globs['tmpdir'] = os.path.join(globs['outdir'], 'tmp');
		else:
			globs['tmpdir'] = args.tmp_dir;
		if not os.path.isdir(globs['tmpdir']):
			os.system("mkdir " + globs['tmpdir']);
	# tmp dir

	globs['sample-name'] = os.path.basename(os.path.normpath(globs['outdir']));
	globs['logfilename'] = os.path.join(globs['outdir'], globs['sample-name'] + ".log");
	if globs['dryrun']:
		globs['logfilename'] = globs['logfilename'].replace(".log", "-dryrun.log");
	globs['scaffs'] = os.path.join(globs['outdir'], "scaffold-list.txt");
	globs['endprog'] = True;
	# Output prep.

	if args.keep_all and args.keep_only_final:
		PC.errorOut("OP8", "Only one of --keepall and --keeponlyfinal can be set.", globs);
	elif args.keep_all:
		globs['keeplevel'] = 2;
	elif args.keep_only_final:
		globs['keeplevel'] = 0;
	# Intermediate file retention options. Keep level 1: Default (keep only final files per iteration). level 2: Keep ALL intermediate files. level 3: Keep only final files from last iteration.

	if args.noindels:
		globs['indels'] = False;
	# Indel output option.

	if args.diploid:
		globs['diploid'] = True;
	# Diploid output option.

	# if args.debug_opt:
	# 	globs['debug'] = True;
	if args.nolog_opt:
		globs['log-v'] = -1;
	# Hidden test options

	globs['num-iters'] = PC.isPosInt(args.num_iters);
	if globs['num-iters'] < 1:
		PC.errorOut("OP9", "-i must be an integer value greater than 1.", globs);
	# Checking the number of iterations option.

	if args.map_only:
		globs['num-iters'] = 1;
		globs['map-only'] = True;
	# Check the map only option.

	if args.filter:
		globs['filter'] = args.filter;
	globs['filter'] = globs['filter'][:-1] + " || ALT=\"*\"'";
	# Check the filter option.

	globs['num-procs'] = PC.isPosInt(args.processes);
	if not globs['num-procs']:
		PC.errorOut("OP10", "-p must be an integer value greater than 1.", globs);
	# Checking the number of processors option.

	if not args.bam:
		if args.map_threads:
			globs['map-t'] = PC.isPosInt(args.map_threads);
			if not globs['map-t']:
				PC.errorOut("OP11", "-map-t must be an integer value greater than 1.", globs);
		elif globs['num-procs'] != 1:
			globs['map-t'] = math.floor(globs['num-procs'] / globs['num-libs']);
	# Getting the number of mapping threads.

	if args.gatk_threads:
		globs['gatk-t'] = PC.isPosInt(args.gatk_threads);
		if not globs['gatk-t']:		
			PC.errorOut("OP12", "-gatk-t must be an integer value greater than 1.", globs);
	if globs['map-only']:
		globs['gatk-t'] = 1;
	# Getting the number of GATK HaplotypeCaller threads

	if globs['num-procs'] > 20:
		globs['gvcf-procs'] = 20;
		globs['filter-procs'] = 20;
	else:
		globs['gvcf-procs'] = globs['num-procs'];
		globs['filter-procs'] = globs['num-procs']
	# Check if the number of proces requested is over the max allowed for GenotypeGVCFs and bcftools filter.

	if globs['num-procs'] < globs['gatk-t'] or globs['num-procs'] < globs['map-t']:
		PC.errorOut("OP13", "-p must be greater than both -bwa-t and -gatk-t, else we can't spawn a single process efficiently.", globs);
	# Check that the procs requested between programs are compatible

	procs_needed = globs['num-libs'] * globs['map-t'];
	if procs_needed > globs['num-procs']:
		globs['map-procs'] = math.floor(globs['num-procs'] / globs['map-t']);
	else:
		globs['map-procs'] = globs['num-libs'];
    # Determine number of mapping processes to launch

	globs['gatk-procs'] = math.floor(globs['num-procs'] / globs['gatk-t']);
    # Determine number of GATK HaplotypeCaller processes to launch

	if args.quiet_flag:
		globs['quiet'] = True;
	# Check the quiet option

	startProg(globs);
	# After all the essential options have been set, call the welcome function.

	if args.quiet_flag:
		globs['log-v'] = 3;
	# If the quiet option was set before, set the verbosity for printWrite here so nothing is
	# printed to the screen, but still output to the logfile.

	return globs;

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

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Output directory:", pad) + globs['outdir']);
	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Temporary file directory:", pad) + globs['tmpdir']);
	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Log file:", pad) + os.path.basename(globs['logfilename']));

	PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
	PC.printWrite(globs['logfilename'], globs['log-v'], "# DEPENDENCY PATHS:");	
	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Program", pad) + "Specified Path");
	# PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
	
	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -mapper", pad) + globs['mapper-path']);
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

	if globs['map-only']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --maponly", pad) +
					PC.spacedOut("True", opt_pad) + 
					"Pseudo-it will do only one iteration and stop after mapping.");

	if globs['bam']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -bam", pad) +
					PC.spacedOut(globs['bam'], opt_pad) + 
					"Pseudo-it will use this BAM file for the first iteration.");

	if globs['keeplevel'] == 0:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --keeponlyfinal", pad) +
					PC.spacedOut("True", opt_pad) + 
					"Pseudo-it will keep only the final files from the last iteration.");
	elif globs['keeplevel'] == 1:			
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --keeponlyfinal, --keepall", pad) +
					PC.spacedOut("False,False", opt_pad) + 
					"Default behavior: Pseudo-it will keep the final files from each iteration.");
	elif globs['keeplevel'] == 2:			
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --keepall", pad) +
					PC.spacedOut("True", opt_pad) + 
					"Pseudo-it will keep ALL files from each iteration.");
	# Reporting the intermediate files options.

	if globs['overwrite']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --overwrite", pad) +
					PC.spacedOut("True", opt_pad) + 
					"Pseudo-it will OVERWRITE the existing files in the specified output directory.");
	# Reporting the overwrite option.

	if globs['resume']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -resume", pad) +
					PC.spacedOut("True", opt_pad) + 
					"Pseudo-it will attempt to resume the run from the specified directory.");
	# Reporting the resume option.

	if not globs['mkdups']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --nomkdup", pad) +
					PC.spacedOut("True", opt_pad) + 
					"Pseudo-it will NOT mark duplicate reads in BAM files with picard MarkDuplicates.");
	else:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --nomkdup", pad) +
					PC.spacedOut("False", opt_pad) + 
					"Pseudo-it will mark duplicate reads in BAM files with picard MarkDuplicates.");		
	# Reporting the mkdups option.	

	if not globs['map-only']:
		if globs['indels']:
			PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --noindels", pad) + 
						PC.spacedOut("False", opt_pad) + 
						"Final assembly will incorporate indels.");
		else:
			PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --noindels", pad) + 
						PC.spacedOut("True", opt_pad) + 
						"Final assembly will NOT incorporate indels.");		
		# Reporting --noindels option.

		if globs['diploid']:
			PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --diploid", pad) + 
						PC.spacedOut("True", opt_pad) + 
						"Final assembly will use IUPAC ambiguity codes for variant sites.");
		else:
			PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --diploid", pad) + 
						PC.spacedOut("False", opt_pad) + 
						"Final assembly will NOT use IUPAC ambiguity codes for variant sites.");		
		# Reporting --diploid option.

		if globs['softmask']:
			PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --hardmask", pad) + 
						PC.spacedOut("False", opt_pad) + 
						"Filtered sites in final consensus will be soft-masked (atcg).");
		else:
			PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --hardmask", pad) + 
						PC.spacedOut("True", opt_pad) + 
						"Filtered sites in final consensus will be hard-masked (Ns).");		
		# Reporting --diploid option.

		if globs['in-vcf']:
			PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -vcf", pad) + 
						PC.spacedOut(globs['in-vcf'], opt_pad) + 
						"Ignoring SNPs in file.");
		# Reporting -vcf option.		

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -i", pad) + 
				PC.spacedOut(str(globs['num-iters']), opt_pad) + 
				"Pseudo-it will perform this many iterations of mapping.");
	# Reporting the number of iterations.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -p", pad) + 
				PC.spacedOut(str(globs['num-procs']), opt_pad) + 
				"Pseudo-it will use this many processes to run.");
	# Reporting the processes option.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -map-t", pad) + 
				PC.spacedOut(str(globs['map-t']), opt_pad) + 
				"pseudo-it will use this many threads for mapping.");
	# Reporting the mapping threads option.

	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Mapping parallel runs", pad) + 
				PC.spacedOut(str(globs['map-procs']), opt_pad) + 
				"This many libraries will be mapped in parallel, each using " + str(globs['map-t']) + " threads.");
	# Reporting the determined number of processes to spawn for mapping.

	if not globs['map-only']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -gatk-t", pad) + 
					PC.spacedOut(str(globs['gatk-t']), opt_pad) + 
					"HaplotypeCaller's --native-pair-hmm-threads option will use this many threads.");
		# Reporting the GATK HaplotypeCaller --native-pair-hmm-threads option.

		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# GATK HaplotypeCaller procs", pad) + 
					PC.spacedOut(str(globs['gatk-procs']), opt_pad) + 
					"This many scaffolds will be called in parallel, each using " + str(globs['gatk-t']) + " threads.");
		# Reporting the determined number of processes to spawn for GATK HaplotypeCaller

		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# GATK GenotypeGVCFs procs", pad) + 
			PC.spacedOut(str(globs['gvcf-procs']), opt_pad) + 
			"This many scaffolds will be genotyped in parallel in the final iteration.");
		# Reporting the determined number of processes to spawn for GATK GenotypeGVCFs

		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# -f", pad) + 
					PC.spacedOut(str(globs['filter']), opt_pad) + 
					"\tVariants will be filtered based on these criteria during the Variant Filtration step.");
		# Reporting the filter option.

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

	# if globs['debug']:
	# 	PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --debug", pad) + 
	# 				PC.spacedOut("True", opt_pad) + 
	# 				"Printing out a bit of debug info.");
	# Reporting the debug option.

	if globs['norun']:
		PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# --norun", pad) + 
					PC.spacedOut("True", opt_pad) + 
					"ONLY PRINTING RUNTIME INFO.");
		PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
	# Reporting the norun option.

#############################################################################