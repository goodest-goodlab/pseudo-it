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
		endProg(globs);
	else:
		sys.exit();

#############################################################################

def fileCheck(globs):
# Checks file options.
	for opt in ['se', 'pe1', 'pe2', 'pem', 'ref']:
		if globs[opt] and not os.path.isfile(globs[opt]):
			errorOut("CORE1", "File not found: " + globs[opt], globs);

#############################################################################

def execCheck(globs, a):
# Checks dependency executables.
	if a.bwa_path:
		globs['bwa-path'] = a.bwa_path;
	if a.picard_path:
		globs['picard-path'] = a.picard_path;
	#globs['picard-path'] = "java -jar " + globs['picard-path'];
	if a.samtools_path:
		globs['samtools-path'] = a.samtools_path;
	if a.gatk_path:
		globs['gatk-path'] = a.gatk_path;
	if a.bedtools_path:
		globs['bedtools-path'] = a.bedtools_path;
	if a.bcftools_path:
		globs['bcftools-path'] = a.bcftools_path;

	for opt in ['bwa-path', 'picard-path', 'samtools-path', 'gatk-path', 'bedtools-path', 'bcftools-path']:
		try:
			cmd_result = subprocess.run(globs[opt], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
		except:
			prog = opt[:opt.index('-')];
			if prog in ['bwa','gatk']:
				prog = prog.upper();
			elif prog == 'picard':
				prog = prog.title();
			errorOut("CORE2", prog + " not found at specified path: " + globs[opt], globs);
	
	return globs;

#############################################################################

def detectRefExt(ref, globs):
# Try and identify the extension of the input FASTA file.
	possible_exts = [".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fna.gz", ".faa", ".faa.gz"];
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
		prev_iter_str = str(globs['iteration'] - 1);
		if len(prev_iter_str) < 2:
			prev_iter_str = "0" + prev_iter_str;
		prev_iter_str = "iter-" + prev_iter_str;
		cur_ref = os.path.join(globs['outdir'], prev_iter_str, prev_iter_str + "-final.fa");
	return cur_ref;

#############################################################################

def runCheck(cur_files, cur_logfile, globs):
# Check whether to run a command or not based on input options and presence of files.

	if globs['overwrite']:
		return True;
	# If overwrite is on, run the command no matter what.

	if globs['resume']:
	# If the output files are present and resume is on:
		if not os.path.isfile(cur_logfile):
		# Check for the log file. If it doesn't exist, re-run the command.
			return True;
		else:
		# Check for the log file. If it does exist, get the last line.
			if os.stat(cur_logfile).st_size == 0:
				return True;
			# If the log file is empty, re-run this step.
			else:
				lastlog = open(cur_logfile, "r").readlines()[-1];
				if "PSEUDOIT SUCCESS!" in lastlog:
					return False;
				# If the Pseudo-it success code is the last line of the file, don't run the command.
				else:
					return True;
				# If the Pseudo-it success code is not the last line of the file, run the command.
	## This version only considers existence of the logfile and the success code. This will be helpful for deleting
	## intermediate files.

	# if all(os.path.isfile(f) for f in cur_files):
	# # If overwrite is not on, action depends on presence of output files.
	# 	if globs['resume']:
	# 	# If the output files are present and resume is on:
	# 		if not os.path.isfile(cur_logfile):
	# 		# Check for the log file. If it doesn't exist, re-run the command.
	# 			return True;
	# 		else:
	# 		# Check for the log file. If it does exist, get the last line.
	# 			lastlog = open(cur_logfile, "r").readlines()[-1];
	# 			if "PSEUDOIT SUCCESS!" in lastlog:
	# 				return False;
	# 			# If the Pseudo-it success code is the last line of the file, don't run the command.
	# 			else:
	# 				return True;
	# 			# If the Pseudo-it success code is not the last line of the file, run the command.
	# 	else:
	# 		return False;
	# 	# This shouldn't happen? If overwrite and resume are off, these files shouldn't exist 
	# 	# (because opt_parse would've stopped the user from using an existing directory).
	## This version considers existence of the actual ouput file.

	else:
		return True;
	# If overwrite isn't on but the files doesn't exist, run the command.

#############################################################################

def runCMD(cmd, cmd_str, cmd_log, report_success, globs):
# Run a command and deal with the output nicely.
	printWrite(globs['logfilename'], globs['log-v'], "# " + getDateTime() + " --> Executing   : " + cmd);
	cmd_result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
	if any(ecode in cmd_result.stderr.decode() for ecode in ['error', 'Error', 'ERROR', 'Exception', 'Could not build fai index', 
																'AssertionError', "Can't read file", "Killed", "No such file or directory", 
																"Symbolic alleles other than <DEL> are currently not supported",
																"Failes to open"]):
		printWrite(globs['logfilename'], globs['log-v'], "\n# " + getDateTime() + " * CMD ERROR: The following command returned an error:\n\n");
		printWrite(globs['logfilename'], globs['log-v'], cmd);
		printWrite(globs['logfilename'], globs['log-v'], "\n\nPlease check the log file for more info: " + cmd_log + "\n");
		printWrite(cmd_log, 3, "CMD:");
		printWrite(cmd_log, 3, cmd + "\n");
		printWrite(cmd_log, 3, "STDOUT output:");
		printWrite(cmd_log, 3, cmd_result.stdout.decode() + "\n");
		printWrite(cmd_log, 3, "STDERR output:");
		printWrite(cmd_log, 3, cmd_result.stderr.decode() + "\n");
		#endProg(globs);
		return True;
	elif report_success:
		printWrite(globs['logfilename'], globs['log-v'], spacedOut("# " + getDateTime() + " --> SUCCESS     : " + cmd_str + " logfile", globs['pad'], sep=".") + cmd_log);
		#printWrite(globs['logfilename'], globs['log-v'], "# " + cmd);
		#printWrite(globs['logfilename'], globs['log-v'], "# Command logfile: " + cmd_log + "\n");
		printWrite(cmd_log, 3, "CMD:");
		printWrite(cmd_log, 3, cmd + "\n");
		printWrite(cmd_log, 3, "STDOUT output:");
		printWrite(cmd_log, 3, cmd_result.stdout.decode() + "\n");
		printWrite(cmd_log, 3, "STDERR output:");
		printWrite(cmd_log, 3, cmd_result.stderr.decode() + "\n");
		printWrite(cmd_log, 3, "PSEUDOIT SUCCESS!");
		return False;
	else:
		printWrite(globs['logfilename'], globs['log-v'], "\n" + getDateTime() + " * 00P5: Something went wrong. This shouldn't have happened. :(\n\n");
		return True;

#############################################################################

def exitCheck(eflag, globs):
	if eflag:
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
	printWrite(globs['logfilename'], globs['log-v'], "# Output directory for this run:   " + globs['outdir']);
	printWrite(globs['logfilename'], globs['log-v'], "# Log file for this run:           " + globs['logfilename']);
	printWrite(globs['logfilename'], globs['log-v'], "# " + "=" * 100);
	printWrite(globs['logfilename'], globs['log-v'], "#");
	sys.exit();

#############################################################################

def getLogTime():
# Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%I.%M.%S");

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
# Add a 0 to a single digit integer....
	iter_str = str(globs['iteration']);
	if len(iter_str) < 2:
		iter_str = "0" + iter_str;
	globs['iter-str'] = iter_str;
	return globs;

#############################################################################

def report_stats(globs, msg="", step_start=0, stat_start=False, stat_end=False, sep=" "):
# Uses psutil to gather memory and time info between steps and print them to the screen.
	import timeit
	if globs['psutil']:
		import psutil;
		dashes = 161;
	else:
		dashes = 101;
	cur_time = timeit.default_timer();
	if stat_start:
	# The first time through just print the headers.
		globs['progstarttime'] = cur_time;
		printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * dashes);
		if globs['psutil']:
			printWrite(globs['logfilename'], globs['log-v'], "# Date/time" + " " * 13 + "Current step" + " " * 25 + "Time since prev (sec)" + " " * 6 + "Elapsed time (sec)" + " " * 4 + "Current mem usage (MB)" + " " * 4 + "Virtual mem usage (MB)");
		else:
			printWrite(globs['logfilename'], globs['log-v'], "# Date/time" + " " * 13 + "Current step" + " " * 25 + "Time since prev (sec)" + " " * 6 + "Elapsed time (sec)");
		printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * dashes);
	else:
		prog_elapsed = round(cur_time - globs['progstarttime'], 5);
		step_elapsed = round(cur_time - step_start, 5);
		if globs['psutil']:
			mem = round(sum([p.memory_info()[0] for p in globs['pids']]) / float(2 ** 20), 5);
			vmem = round(sum([p.memory_info()[1] for p in globs['pids']]) / float(2 ** 20), 5);
			printWrite(globs['logfilename'], globs['log-v'], "# " + getDateTime() + " " + msg + sep * (37-len(msg)) + str(step_elapsed) + sep * (27-len(str(step_elapsed))) + str(prog_elapsed) + sep * (22-len(str(prog_elapsed))) + str(mem) + sep * (26-len(str(mem))) + str(vmem));
		else:
			printWrite(globs['logfilename'], globs['log-v'], "# " + getDateTime() + " " + msg + sep * (37-len(msg)) + str(step_elapsed) + sep * (27-len(str(step_elapsed))) + str(prog_elapsed) + sep * (22-len(str(prog_elapsed))));
		if stat_end:
			printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * dashes);
	return cur_time;

#############################################################################

def getSubPID(n):
# Gets the process ids for the --stats option.
	import psutil
	return psutil.Process(os.getpid());

#############################################################################

def welcome():
# Reads the ASCII art "Referee" text to be printed to the command line.
	return open(os.path.join(os.path.dirname(__file__), "pi-welcome.txt"), "r").read();

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
