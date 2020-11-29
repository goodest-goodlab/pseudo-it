#############################################################################
# Pseudo-it CORE functions
# Gregg Thomas
# August 2013-present
# Forked from Referee on 12.07.2019
# Updated for Pseudo-it December 2019
#############################################################################

import sys, os, timeit, datetime, time, gzip, string, random, subprocess

#############################################################################

def errorOut(errnum, errmsg, globs):
# Formatting for error messages.
	fullmsg = "**Error " + str(errnum) + ": " + errmsg;
	border = "-" * len(fullmsg);
	fullstr = "\n" + border + "\n" + fullmsg + "\n" + border + "\n"
	printWrite(globs['logfilename'], globs['log-v'], "\n" + border + "\n" + fullmsg + "\n" + border + "\n");
	if globs['endprog']:
		globs['exit-code'] = 1;
		endProg(globs);
	else:
		printWrite(globs['logfilename'], globs['log-v'], "\nScript call: " + " ".join(sys.argv));
		sys.exit(1);

#############################################################################

def fileCheck(globs):
# Checks file options.
	files = ['se', 'pe1', 'pe2', 'pem', 'ref'];
	if globs['bam']:
		files += ['bam', 'bam-index'];
	for opt in files:
		if globs[opt] and not os.path.isfile(globs[opt]):
			errorOut("CORE1", "File not found: " + globs[opt], globs);

#############################################################################

def execCheck(globs, a):
# Checks dependency executables.
	deps_passed = True;
	# Variable to check if all dependencies are found.

	if a.bwa_path:
		globs['bwa-path'] = a.bwa_path;
	if a.picard_path:
		globs['picard-path'] = a.picard_path;
	if a.samtools_path:
		globs['samtools-path'] = a.samtools_path;
	if a.gatk_path:
		globs['gatk-path'] = a.gatk_path;
	if a.bedtools_path:
		globs['bedtools-path'] = a.bedtools_path;
	if a.bcftools_path:
		globs['bcftools-path'] = a.bcftools_path;
	# Update the global paths if the user provided them through args.

	dpad = 14;
	if a.depcheck:
		print("# --depcheck set: CHECKING DEPENDENCY PATHS AND EXITING.\n");
		print(spacedOut("   PROGRAM", dpad) + spacedOut("PATH", dpad) + "STATUS");
		print("   -------------------------------");
	# For the dependency check option (--depcheck), this initializes a neat output table.

	for opt in ['bwa-path', 'picard-path', 'samtools-path', 'gatk-path', 'bedtools-path', 'bcftools-path', 'tabix-path']:
		
		prog = opt[:opt.index('-')];
		if prog in ['bwa','gatk']:
			prog = prog.upper();
		elif prog == 'picard':
			prog = prog.title();
		# Get a formatted program name.

		dcheck_str = [spacedOut("   " + prog, dpad), spacedOut(globs[opt], dpad), "NA"];
		# Initialize the check string for --depcheck.

		cmd_result = subprocess.run(globs[opt], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
		# Run the provided command and retrieve the exit code.

		if cmd_result.returncode > 1:
		# If the exit code for the command run is greater than 1, the command isn't found.
			dcheck_str[2] = "FAILED with exit code " + str(cmd_result.returncode);
			deps_passed = False;
			# Update the check string and keep going.
			if not a.depcheck:	
				errorOut("CORE2", prog + " not found at specified path: " + globs[opt], globs);
			# On a normal run, exit immediately.
		else:
			dcheck_str[2] = "PASSED";
			# Update the check string.
			
		if a.depcheck:
			print("".join(dcheck_str));
		# Print the check string if --depcheck is set.

	return globs, deps_passed;

#############################################################################

def detectRefExt(ref, globs):
# Try and identify the extension of the input FASTA file.
	possible_exts = [".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fna.gz"];
	for e in possible_exts:
		if ref.endswith(e):
			return e;
	errorOut("CORE3", "Make sure your reference FASTA (-ref) has one of the following extensions: " + ", ".join(possible_exts), globs);

#############################################################################

def getRef(globs):
# Get the reference file based on the iteration number.
	if globs['iteration'] == 1:
		cur_ref = globs['ref'];
	else:
		cur_ref = globs['consensus-file'];
	return cur_ref;

#############################################################################

def runCheck(cmd, cmds, globs):
# Check whether to run a command or not based on input options and presence of files.
	if os.path.isfile(cmds[cmd]['outfile']) and os.stat(cmds[cmd]['outfile']).st_size != 0:
		if os.path.isfile(cmds[cmd]['logfile']) and os.stat(cmds[cmd]['logfile']).st_size != 0:
			log_last_line = open(cmds[cmd]['logfile'], "r").readlines()[-1];
			if "PSEUDOIT SUCCESS!" in log_last_line:
				report_step(globs, cmds, cmd, "RESUME", "previous output found: " + cmds[cmd]['logfile']);
				return False;
	# If the output file exists and has content AND the logfile exists and has content AND the last line of the logfile is the pseudo-it success code, 
	# then return False to NOT run the command
	return True;
	# In every other case, return True to run the command

#############################################################################

def prevCheck(outfile, logfile, globs):
	if os.path.isfile(outfile) and os.stat(outfile).st_size != 0:
		if os.path.isfile(logfile) and os.stat(logfile).st_size != 0:
			log_last_line = open(logfile, "r").readlines()[-1];
			if "PSEUDOIT SUCCESS!" in log_last_line:
				return False;
	# If the output file exists and has content AND the logfile exists and has content AND the last line of the logfile is the pseudo-it success code, 
	# then return False to NOT run the command
	return True;
	# In every other case, return True to run the command	

#############################################################################

def runCMD(cmd, globs, cmds, report_success):
# Run a command and deal with the output nicely
	if globs['dryrun']:
		report_step(globs, cmds, cmd, "DRYRUN", cmd);
		return False;
	# If this is a dry run, don't run the command.

	ecodes = ['error', 'Error', 'ERROR', 'Exception', 'Could not build fai index', 
				'AssertionError', "Can't read file", "Killed", "No such file or directory", 
				"Symbolic alleles other than <DEL> are currently not supported",
				"Failed to open", "The index file is older than the data file",
				"command not found"]

	run = True;
	if globs['resume']:
		run = runCheck(cmd, cmds, globs);
	# If -resume is set, check some conditions to determine whether to run the command

	if run:
		report_step(globs, cmds, cmd, "EXECUTING", cmd);
		cmd_result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
		cmd_stdout = cmd_result.stdout.decode();
		cmd_stderr = cmd_result.stderr.decode();
		cmd_output = cmd_stdout + "\n\n" + cmd_stderr;

		if "The index file is older than the data file" in cmd_output:
			reindex_cmd = "tabix -fp vcf " + cmds[cmd]['vcffile'];
			report_step(globs, cmds, cmd, "Re-indexing", reindex_cmd);
			cmd_result = subprocess.run(reindex_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
			if any(ecode in cmd_result.stderr.decode() for ecode in ecodes):
				report_step(globs, cmds, cmd, "RE-INDEX ERROR!", "Check log: " + cmds[cmd]['logfile']);
				return True;
			else:
				re_run_status = runCMD(cmd, globs, cmds, report_success);

			if re_run_status:
				return True;
			else:
				return False;
		# This block handles re-indexing for VCF files if GATK indexes it too quickly.

		printWrite(cmds[cmd]['logfile'], 3, "CMD:");
		printWrite(cmds[cmd]['logfile'], 3, cmd + "\n");
		printWrite(cmds[cmd]['logfile'], 3, "STDOUT output:");
		printWrite(cmds[cmd]['logfile'], 3, cmd_stdout + "\n");
		printWrite(cmds[cmd]['logfile'], 3, "STDERR output:");
		printWrite(cmds[cmd]['logfile'], 3, cmd_stderr + "\n");

		if any(ecode in cmd_output for ecode in ecodes):
			report_step(globs, cmds, cmd, "ERROR!", "Check log: " + cmds[cmd]['logfile']);
			return True;
		elif report_success:
			printWrite(cmds[cmd]['logfile'], 3, "PSEUDOIT SUCCESS!");
			report_step(globs, cmds, cmd, "SUCCESS", cmds[cmd]['logfile']);
			return False;		
		else:
			printWrite(globs['logfilename'], globs['log-v'], "\n" + getDateTime() + " * 00P5: Something went wrong. This shouldn't have happened. :(\n\n");
			return True;
	else:
		return False;

#############################################################################

def getLogTime():
# Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%I.%M.%S");

#############################################################################

def getDate():
# Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y");

#############################################################################

def getTime():
# Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%H:%M:%S");

#############################################################################

def getDateTime():
# Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y | %H:%M:%S");

#############################################################################

def getOutTime():
# Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m-%d-%Y.%I-%M-%S");

#############################################################################

def printWrite(o_name, v, o_line1, o_line2="", pad=0):
# Function to print a string AND write it to the file.
	if o_line2 == "":
		outline = o_line1;
	else:
		outline = o_line1 + " "*(pad-len(o_line1)) + o_line2;
	if v in [-1,1,2]:
		print(outline);
	if v != -1:
		f = open(o_name, "a");
		f.write(outline + "\n");
		f.close();

#############################################################################
	
def spacedOut(string, totlen, sep=" "):
# Properly adds spaces to the end of a message to make it a given length
	spaces = sep * (totlen - len(string));
	return string + spaces;

#############################################################################

def getIterStr(globs):
# Add a 0 to a single digit integer so they align with double digit integers when printed
	iter_str = str(globs['iteration']);
	if len(iter_str) < 2:
		iter_str = "0" + iter_str;
	globs['iter-str'] = iter_str;
	return globs;

#############################################################################

def getCMDNum(globs, num_cmds):
# Format the command numbers for the log
	cmd_num = str(num_cmds+1);
	while len(cmd_num) < 4:
		cmd_num = "0" + cmd_num;
	cmd_num = globs['iter-str'] + "-" + cmd_num;
	return cmd_num;

#############################################################################

def isPosInt(numstr):
# Check if a string is a positive integer
	try:
		num = int(numstr);
	except:
		return False;

	if num > 0:
		return num;
	else:
		return False;

#############################################################################

def report_step(globs, cmds, cmd="", status="", result="", start=False, end=False, sep=" "):
	cur_time = timeit.default_timer();
	dashes = 150;
	col_widths = [ 14, 10, 20, 20, 50, 12, 25 ];
	if start:
		headers = [ "# Date", "Time", "Elapsed time (s)", "Step time (s)", "Current step", "Status"];
		headers = "".join([ spacedOut(str(headers[i]), col_widths[i]) for i in range(len(headers)) ]);

		print("# " + "-" * 125);
		print(headers);
		print("# " + "-" * 125);
		
		printWrite(globs['logfilename'], 3, "# " + "-" * dashes);
		headers += "Command/Result";
		printWrite(globs['logfilename'], 3, headers);
		printWrite(globs['logfilename'], 3, "# " + "-" * dashes);

		return cur_time;

	elif end:
		iter_elapsed = str(round(cur_time - globs['iterstarttime'], 5));
		#prog_elapsed = str(round(cur_time - globs['progstarttime'], 5));

		#printWrite(globs['logfilename'], globs['log-v'], "# Program elapsed time (s):   " + prog_elapsed);
		printWrite(globs['logfilename'], globs['log-v'], "# Iteration elapsed time (s): " + iter_elapsed);

	else:
		prog_elapsed = str(round(cur_time - globs['progstarttime'], 5));
		if cmd.startswith("NA"):
			step_elapsed = "-";
			step = globs['iter-str'] + cmd.replace("NA-", "");
		else:
			if cmds[cmd]['start']:
				step_elapsed = str(round(cur_time - cmds[cmd]['start'], 5));
			else:
				cmds[cmd]['start'] = cur_time;
				step_elapsed = "-";
			step = cmds[cmd]['cmd-num'] + " " + cmds[cmd]['desc'];
		printline = [ "# " + getDate(), getTime(), prog_elapsed, step_elapsed, step, status ];
		printline = [ spacedOut(str(printline[i]), col_widths[i]) for i in range(len(printline)) ];
		print("".join(printline));

		logline = [ "# " + getDate(), getTime(), prog_elapsed, step_elapsed, step, status, result ];
		logline = [ spacedOut(str(logline[i]), col_widths[i]) for i in range(len(logline)) ];
		printWrite(globs['logfilename'], 3, "".join(logline));

	#if end:
		#printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * dashes);


#############################################################################

def welcome():
# Reads the ASCII art "Referee" text to be printed to the command line.
	return open(os.path.join(os.path.dirname(__file__), "pi-welcome.txt"), "r").read();

#############################################################################

def exitCheck(eflag, globs):
	if eflag:
		globs['exit-code'] = 1;
		endProg(globs);

#############################################################################

def endProg(globs):
# A nice way to end the program.
	if globs['quiet']:
		globs['log-v'] = 1;
	endtime = timeit.default_timer();
	totaltime = endtime - globs['starttime'];

	printWrite(globs['logfilename'], globs['log-v'], "#\n# Done!");
	printWrite(globs['logfilename'], globs['log-v'], "# The date and time at the end is: " + getDateTime());
	printWrite(globs['logfilename'], globs['log-v'], "# Total execution time:            " + str(round(totaltime,3)) + " seconds.");
	if globs['exit-code'] == 0:
		if globs['map-only']:
			printWrite(globs['logfilename'], globs['log-v'], "# Final BAM file:                  " + globs['iter-final-bam']); 
		else:
			printWrite(globs['logfilename'], globs['log-v'], "# Final Assembly:                  " + globs['consensus-file']);
	
	printWrite(globs['logfilename'], globs['log-v'], "# Output directory for this run:   " + globs['outdir']);
	printWrite(globs['logfilename'], globs['log-v'], "# Log file for this run:           " + globs['logfilename']);

	if globs['exit-code'] != 0:
		printWrite(globs['logfilename'], globs['log-v'], "#\n# ERROR: NON-ZERO EXIT STATUS.");
		printWrite(globs['logfilename'], globs['log-v'], "# ERROR: PSEUDO-IT FINISHED WITH ERRORS.");
		printWrite(globs['logfilename'], globs['log-v'], "# ERROR: PLEASE CHECK THE LOG FILE FOR MORE INFO: " + globs['logfilename'] + "\n#");

	print("# " + "=" * 125);
	printWrite(globs['logfilename'], 3, "# " + "=" * 150);
	printWrite(globs['logfilename'], globs['log-v'], "#");
	sys.exit(globs['exit-code']);

#############################################################################














































































































































































































































































































































































































































































































































































































































































































































































def simpson():
	s = """
		              @                                                           
	             CC   CQ                                                      
	            /CCB @CC                                                      
	        GCCS CCCCCCC7                                                     
	         @CCCCCCCC@@@                                                     
	        @@@CCCCCCCCCCCCC                                                  
	        @CCCCCCCCCCCCCCCCC/                                               
	          OCCCCCCCCCCCCCCCCC                                              
	          @CCCCCCCCCCCCCCCCCCC@                                           
	          CCCCCCCCCCC@QCCCCCC@CCC(                                        
	          6CCCCCCCCCCCOCCCCCCCCC@es@                                      
	          @CCCCCCCCCCCCCCCCCCCCCCCCCC                                     
	          ^CCCCCCCCCCCCC@      @K      R                                  
	           CCCC@CCCCCCC                                                   
	          ~CCCC@CCCCCC@          G     //                   SCC@  @CC~    
	          @CCCCsCCCCCC#    S@    #       RS@                CCCB 6CCC     
	          @CCCCCCCCCOOG        S/@OC@CCSR  @/K             GCCCC@CCCC @CC 
	           @CC@GCCCCCCCCS      K67@CCCCCCCCG @             CCCCKCCCC@CCCC 
	          sCCCCCCCCCCCCCCCCCC77777SCCCCCCCCCCC        @@   CCCCKCCCCCCCCK 
	          CCCGRCCCCCCCCCCCCCCC777@CC@CCCCCCCCC       @CCCC@CCCCCCCCCCCCR  
	          @CCCK@CCCCCCCCCCCCCCCSQ(((((((((((@         6CCCCCCCCCCCCCCCC(  
	           ^QCCCCCCCCCCCCCCCC%(((((((((((((((((%@s     @CCCCCCCCCCCCCCC   
	            3CBCCCCCCCCCCR(((((((((((((((((((((((((%    CCCCCCC@CCCCCCC   
	           /CCCCCCCCCCC(((((((((((((((((((@((((@((@     SCCCCCCCCCCCCCC   
	           6CCCC@@CCCC(((((@   #@@@@       @@K          @KCCCCCC@CCCCCC   
	           CCCCCBCCCCC(((((@@@@@@@@(((3                @((#CCCCCKCCCCC    
	           CCCCSCCCCCC@((((@KK@KR@(((               %~(%(((%sCCCCCCCC     
	           ~@CCCCCCC@CCC((((((O@@@%((C             @~(((@(((((~(6((%      
	          (((((@CCCCRCCC@(((((((((((@@            @(((((((@t(((((((       
	          ((((((((@@@@CCC#(((((((((@%@(@         s(((((((((((((Ct(        
	         %(@(((((((((((((((((((@   /((((       @(((((((((((((((((@        
	        sR(((@((((((((((((((((@      (#%@G((((((((((((((((((((((@         
	      7~((((((((t@((((((((((G(s      ((((((((((((((((((((((((((R          
	     B((((((((((((((s#((((((((((/  s(e((((((((((((((((((((((((R           
	    @((((((((((((((((((((((((((((@((((@((((((((((((((((((((((/            
	   /(((((((((((((((((((((((((((((((((((((((((((((((((((((((C              
	   ~(((((((((((((((((((((((((((((e((((((((((((((((((((((((@               
	  @((((((((((((((((((((((((((((((@(((@(#((((((((((((((((K                 
	  G((((((((((((((((((((((((((((((s((((%@(((((((((((((((6                  
	  (((((((((((((((((((((((((((((((((((((e(((((((((((Q/                     
	 G((((((((((((((@((((((((((((((((((((((%(((((((@                          
	 @(((((((((((((C(((((((((((((((((((t(((G(                                 
	 ~(((((((((((((((((((((((((((((((((t(((e(s                                
	@(((((((((((((((((((((((((((((((((((((((((                                
	7(((((((((((((6((((((((((((((((((((Q((((((%   
	"""

	print(s);
