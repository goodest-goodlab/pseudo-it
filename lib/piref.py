# Reference FASTA functions for Pseudo-it
#############################################################################
import sys, os, multiprocessing as mp, subprocess, lib.picore as PC

#############################################################################

def indexCheck(cur_fa, globs, cmds):
# Checks that the user has created the proper index files before running the program.

    ref_ext = PC.detectRefExt(cur_fa, globs);

    dictfile = cur_fa.replace(ref_ext, ".dict");
    cmd = "os.path.isfile(" + dictfile + ")";
    cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Checking ref indices (.dict)", 'outfile' : "", 'logfile' : "", 'start' : False };
    PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
    if not os.path.isfile(dictfile):
        PC.errorOut("REF1", "Reference dictionary not found. Please run: picard CreateSequenceDictionary R=<ref>.fa O=<ref>.dict", globs);
    PC.report_step(globs, cmds, cmd, "SUCCESS", "index file found (.dict)", "");
    # Check for the reference dictionary file.

    faidxfile = cur_fa + ".fai";
    cmd = "os.path.isfile(" + faidxfile + ")";
    cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Checking ref indices (.fai)", 'outfile' : "", 'logfile' : "", 'start' : False };
    PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
    if not os.path.isfile(faidxfile):
        PC.errorOut("REF2", "Reference index (samtools) not found. Please run: samtools faidx <ref>.fa", globs);
    PC.report_step(globs, cmds, cmd, "SUCCESS", "index file found (.fai)");
    # Check for the reference faidx file.

    if globs['mapper'] == "bwa":    
        indexfiles = [cur_fa + ".amb", cur_fa + ".ann", cur_fa + ".bwt", cur_fa + ".pac", cur_fa + ".sa"];
        cmd = "os.path.isfile(" + ",".join(indexfiles) + ")";
        cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Checking ref indices (" + globs['mapper'] + ")", 'outfile' : "", 'logfile' : "", 'start' : False };
        PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
        if any(not os.path.isfile(f) for f in indexfiles):
            PC.errorOut("REF3", "Reference index (bwa) not found. Please run: bwa index <ref>.fa", globs);
        PC.report_step(globs, cmds, cmd, "SUCCESS", "index files found (" + globs['mapper'] + ")");
    # Check for the bwa index files if --mapper is bwa.

    elif globs['mapper'] == "hisat2":
        indexfile = cur_fa + ".1.ht2";
        cmd = "os.path.isfile(" + indexfile + ")";
        cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Checking ref indices (" + globs['mapper'] + ")", 'outfile' : "", 'logfile' : "", 'start' : False };
        PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
        if not os.path.isfile(indexfile):
            PC.errorOut("REF3", "Reference index (hisat2) not found. Please run: hisat2-build <ref>.fa <ref>.fa", globs);
        PC.report_step(globs, cmds, cmd, "SUCCESS", "index file found (" + globs['mapper'] + ")");
    # Check for the hisat2 index files if --mapper is hisat2.        

    return cmds;

#############################################################################

def getScaffs(cur_fa, globs, cmds, report_status=True):
# Save the list of scaffolds/contigs/chromosomes from a FASTA file to a text file.

    cmd = "getScaffs()";
    cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Get ref scaffold IDs from .fai index", 'outfile' : "", 'logfile' : "", 'start' : False };
    # Add the command to the global commands dict.

    if not globs['dryrun']:
        PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
        indexfile = cur_fa + ".fai";
        scaffs = [ line.strip().split("\t")[0] for line in open(indexfile) ];
        PC.report_step(globs, cmds, cmd, "SUCCESS", str(len(scaffs)) + " scaffold IDs read");
    else:
        PC.report_step(globs, cmds, cmd, "DRYRUN", cmd);
        scaffs = [];       

    return scaffs, cmds;

    # cmd = "grep \">\" " + cur_fa + " | sed 's/>//g'"# > " + globs['scaffs'];
    # # grep the number of scaffolds in the reference... I guess this could also be done by just reading
    # # the number of lines in the index file...

    # cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Get ref scaffold IDs", 'outfile' : "", 'logfile' : "", 'start' : False };
    # # Add the grep command to the global commands dict.

    # if not globs['dryrun']:
    #     PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
    #     cmd_result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
    #     cur_scaffs = list(filter(None, cmd_result.stdout.decode().split("\n")));
    #     globs['scaffolds'] = [ scaff[:scaff.index(" ")] if " " in scaff else scaff for scaff in cur_scaffs ];
    #     PC.report_step(globs, cmds, cmd, "SUCCESS", str(len(globs['scaffolds'])) + " scaffold IDs read");
    # else:
    #     PC.report_step(globs, cmds, cmd, "DRYRUN", cmd);
    #     globs['scaffolds'] = [];
    # # Run the grep command and check for errors..
        
    # return cmds;

#############################################################################

def indexFa(globs, cmds, cur_ref):
# Creates all reference fasta index files for subsequent iterations. For the first
# iteration these are assumed to be created before the program is run.

    indices = ['dict', 'faidx', 'index'];
    # The types of indices needed: .dict from picard, .fai from samtools, and the current --mapper index files.

    index_cmds = {};
    ref_ext = PC.detectRefExt(cur_ref, globs);
    # Detect whether the reference is compressed or not.

    for step in indices:
        if step == 'dict':
            cur_logfile = os.path.join(globs['iterlogdir'], "picard-dict-iter-" + globs['iter-str'] + ".log");
            dict_file = cur_ref.replace(ref_ext, ".dict");

            if os.path.isfile(dict_file) and globs['overwrite']:
                os.system("rm " + dict_file);

            picard_cmd = globs['picard-path'] + " CreateSequenceDictionary R=" + cur_ref + " O=" + dict_file;
            cmds[picard_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create reference dict", 'outfile' : dict_file, 'logfile' : cur_logfile, 'start' : False };
            index_cmds[picard_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create reference dict", 'outfile' : dict_file, 'logfile' : cur_logfile, 'start' : False };
        # Create the reference dictionary by running picard CreateSequenceDictionary

        if step == "faidx":
            cur_logfile = os.path.join(globs['iterlogdir'], "samtools-faidx-iter-" + globs['iter-str'] + ".log");
            faidx_file = cur_ref + ".fai";            

            faidx_cmd = globs['samtools-path'] + " faidx " + cur_ref;
            cmds[faidx_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create reference faidx", 'outfile' : faidx_file, 'logfile' : cur_logfile, 'start' : False };
            index_cmds[faidx_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create reference faidx", 'outfile' : faidx_file, 'logfile' : cur_logfile, 'start' : False };
        # Create the reference index by running samtools faidx

        if step == "index":
            if globs['mapper'] == "bwa":
                cur_logfile = os.path.join(globs['iterlogdir'], "bwa-index-iter-" + globs['iter-str'] + ".log");
                index_files = [cur_ref + ".amb", cur_ref + ".ann", cur_ref + ".bwt", cur_ref + ".pac", cur_ref + ".sa"];

                index_cmd = globs['mapper-path'] + " index " + cur_ref;
                cmds[index_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create BWA reference index", 'outfile' : "", 'logfile' : cur_logfile, 'start' : False };
                index_cmds[index_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create BWA reference index", 'outfile' : "", 'logfile' : cur_logfile, 'start' : False };
            # Create the reference index by running bwa index if --mapper is bwa

            elif globs['mapper'] == "hisat2":
                cur_logfile = os.path.join(globs['iterlogdir'], "hisat2-build-index-iter-" + globs['iter-str'] + ".log");
                index_file = cur_ref + ".ht";

                index_cmd = globs['mapper-path'] + "-build " + cur_ref + " " + cur_ref;
                cmds[index_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create hisat2 reference index", 'outfile' : "", 'logfile' : cur_logfile, 'start' : False };
                index_cmds[index_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create hisat2 reference index", 'outfile' : "", 'logfile' : cur_logfile, 'start' : False };                
            # Create the reference index by running hisat2-build if --mapper is hisat2

    index_procs = min(3, globs['num-procs']);
    pool = mp.Pool(processes=index_procs);
    for result in pool.starmap(PC.runCMD, ((index_cmd, globs, cmds, True) for index_cmd in index_cmds )):
        if result:
            pool.terminate();
            globs['exit-code'] = 1;
            PC.endProg(globs);
    pool.terminate();
    # Run the index commands in parallel and check for errors.

    return cmds;

#############################################################################