# Reference FASTA functions for Pseudo-it
#############################################################################
import sys, os, subprocess, lib.picore as PC
#############################################################################

def getScaffs(cur_fa, globs, step_start_time, report_status=True):
# Save the list of scaffolds/contigs/chromosomes from a FASTA file to a text file.
    if report_status:
        step_start_time = PC.report_stats(globs, "PREP: Get scaffold IDs", step_start=step_start_time);
    cmd = "grep  \">\" " + cur_fa + " | sed 's/>//g'"# > " + globs['scaffs']; 
    cmd_result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
    cur_scaffs = list(filter(None, cmd_result.stdout.decode().split("\n")));
    cur_scaffs = [ scaff[:scaff.index(" ")] if " " in scaff else scaff for scaff in cur_scaffs ];
    return cur_scaffs, step_start_time;

#############################################################################

def indexCheck(cur_fa, globs):
# Checks that the user has created the proper index files before running the program.
    ref_ext = PC.detectRefExt(cur_fa, globs);
    dictfile = cur_fa.replace(ref_ext, ".dict");
    if not os.path.isfile(dictfile):
        PC.errorOut("REF1", "Reference dictionary not found. Please run: picard CreateSequenceDictionary R=<ref>.fa O=<ref>.dict", globs);

    faidxfile = cur_fa + ".fai";
    if not os.path.isfile(faidxfile):
        PC.errorOut("REF2", "Reference index (samtools) not found. Please run: samtools faidx <ref>.fa", globs);
        
    indexfiles = [cur_fa + ".amb", cur_fa + ".ann", cur_fa + ".bwt", cur_fa + ".pac", cur_fa + ".sa"];
    if any(not os.path.isfile(f) for f in indexfiles):
        PC.errorOut("REF3", "Reference index (bwa) not found. Please run: bwa index <ref>.fa", globs);

#############################################################################

def indexDistributor(index_item):
# This function distributes the FASTA indexing jobs to be run in parallel.
    cur_opt, cur_ref, globs, step_start_time = index_item;
    if cur_opt == 'dict':
        step_start_time, exit_flag = picardDict(cur_ref, globs, step_start_time);
    elif cur_opt == 'faidx':
        step_start_time, exit_flag = faidxFasta(cur_ref, globs, step_start_time);
    elif cur_opt == 'index':
        step_start_time, exit_flag = indexFasta(cur_ref, globs, step_start_time);
    
    return step_start_time, exit_flag;

############################################################################

def picardDict(cur_ref, globs, step_start_time):
# Create the sequence dictionary for the FASTA.    
    ref_ext = PC.detectRefExt(cur_ref, globs);
    cur_logfile = os.path.join(globs['logdir'], "picard-dict-iter-" + globs['iter-str'] + ".log");
    dictfile = cur_ref.replace(ref_ext, ".dict");

    run_flag = PC.runCheck([dictfile], cur_logfile, globs);

    if run_flag:
        #step_start_time = PC.report_stats(globs, "--> ITERATION " + globs['iter-str'] + ": Dictionary", step_start=step_start_time);
        dict_cmd = globs['picard-path'] + " CreateSequenceDictionary R=" + cur_ref + " O=" + dictfile;
        exit_flag = PC.runCMD(dict_cmd, "Picard CreateSequenceDictionary", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Dictionary file already exists", globs['pad'], sep=".") + dictfile);
        exit_flag = False;
    return step_start_time, exit_flag;

#############################################################################

def faidxFasta(cur_ref, globs, step_start_time):
# Index the FASTA with samtools faidx.
    cur_logfile = os.path.join(globs['logdir'], "samtools-faidx-iter-" + globs['iter-str'] + ".log");
    faidxfile = cur_ref + ".fai";

    run_flag = PC.runCheck([faidxfile], cur_logfile, globs);

    if run_flag:        
        #step_start_time = PC.report_stats(globs, "--> ITERATION " + globs['iter-str'] + ": faidx", step_start=step_start_time);
        faidx_cmd = globs['samtools-path'] + " faidx " + cur_ref;
        exit_flag = PC.runCMD(faidx_cmd, "samtools faidx", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# faidx file already exists", globs['pad'], sep=".") + faidxfile);
        exit_flag = False;
    # samtools faidx
    return step_start_time, exit_flag;

#############################################################################

def indexFasta(cur_ref, globs, step_start_time):
# Index the FASTA file with BWA.
    cur_logfile = os.path.join(globs['logdir'], "bwa-index-iter-" + globs['iter-str'] + ".log");
    indexfiles = [cur_ref + ".amb", cur_ref + ".ann", cur_ref + ".bwt", cur_ref + ".pac", cur_ref + ".sa"];

    run_flag = PC.runCheck(indexfiles, cur_logfile, globs);

    if run_flag:
        #step_start_time = PC.report_stats(globs, "--> ITERATION " + globs['iter-str'] + ": BWA index", step_start=step_start_time);
        index_cmd = globs['bwa-path'] + " index " + cur_ref;
        exit_flag = PC.runCMD(index_cmd, "BWA index", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# BWA index files already exist", globs['pad'], sep=".") + ",".join(indexfiles));
        exit_flag = False;
    # bwa index
    return step_start_time, exit_flag;

#############################################################################

