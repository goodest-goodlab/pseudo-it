# Reference FASTA functions for Pseudo-it
#############################################################################
import sys, os, multiprocessing as mp, subprocess, lib.picore as PC

#############################################################################

def indexCheck(cur_fa, globs, cmds):
# Checks that the user has created the proper index files before running the program.

    ref_ext = PC.detectRefExt(cur_fa, globs);

    dictfile = cur_fa.replace(ref_ext, ".dict");
    cmd = "os.path.isfile(" + dictfile + ")";
    cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Checking ref indices", 'outfile' : "", 'logfile' : "", 'start' : False };
    PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
    if not os.path.isfile(dictfile):
        PC.errorOut("REF1", "Reference dictionary not found. Please run: picard CreateSequenceDictionary R=<ref>.fa O=<ref>.dict", globs);
    PC.report_step(globs, cmds, cmd, "SUCCESS", "index file found", "");

    faidxfile = cur_fa + ".fai";
    cmd = "os.path.isfile(" + faidxfile + ")";
    cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Checking ref indices", 'outfile' : "", 'logfile' : "", 'start' : False };
    PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
    if not os.path.isfile(faidxfile):
        PC.errorOut("REF2", "Reference index (samtools) not found. Please run: samtools faidx <ref>.fa", globs);
    PC.report_step(globs, cmds, cmd, "SUCCESS", "index file found");
        
    indexfiles = [cur_fa + ".amb", cur_fa + ".ann", cur_fa + ".bwt", cur_fa + ".pac", cur_fa + ".sa"];
    cmd = "os.path.isfile(" + ",".join(indexfiles) + ")";
    cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Checking ref indices", 'outfile' : "", 'logfile' : "", 'start' : False };
    PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
    if any(not os.path.isfile(f) for f in indexfiles):
        PC.errorOut("REF3", "Reference index (bwa) not found. Please run: bwa index <ref>.fa", globs);
    PC.report_step(globs, cmds, cmd, "SUCCESS", "index files found");

    return cmds;

#############################################################################

def getScaffs(cur_fa, globs, cmds, report_status=True):
# Save the list of scaffolds/contigs/chromosomes from a FASTA file to a text file.

    cmd = "grep \">\" " + cur_fa + " | sed 's/>//g'"# > " + globs['scaffs'];
    cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Get ref scaffold IDs", 'outfile' : "", 'logfile' : "", 'start' : False };

    if not globs['dryrun']:
        PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
        cmd_result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
        cur_scaffs = list(filter(None, cmd_result.stdout.decode().split("\n")));
        globs['scaffolds'] = [ scaff[:scaff.index(" ")] if " " in scaff else scaff for scaff in cur_scaffs ];
        PC.report_step(globs, cmds, cmd, "SUCCESS", str(len(globs['scaffolds'])) + " scaffold IDs read");

        return cmds;
    else:
        PC.report_step(globs, cmds, cmd, "DRYRUN", cmd);
        globs['scaffolds'] = [];
        return cmds

#############################################################################

def indexFa(globs, cmds, cur_ref):
# Creates all reference fasta index files

    indices = ['dict', 'faidx', 'index'];
    index_cmds = {};
    ref_ext = PC.detectRefExt(cur_ref, globs);
    for step in indices:
        if step == 'dict':
            cur_logfile = os.path.join(globs['iterlogdir'], "picard-dict-iter-" + globs['iter-str'] + ".log");
            dict_file = cur_ref.replace(ref_ext, ".dict");

            picard_cmd = globs['picard-path'] + " CreateSequenceDictionary R=" + cur_ref + " O=" + dict_file;
            cmds[picard_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create reference dict", 'outfile' : dict_file, 'logfile' : cur_logfile, 'start' : False };
            index_cmds[picard_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create reference dict", 'outfile' : dict_file, 'logfile' : cur_logfile, 'start' : False };

        if step == "faidx":
            cur_logfile = os.path.join(globs['iterlogdir'], "samtools-faidx-iter-" + globs['iter-str'] + ".log");
            faidx_file = cur_ref + ".fai";            

            faidx_cmd = globs['samtools-path'] + " faidx " + cur_ref;
            cmds[faidx_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create reference faidx", 'outfile' : faidx_file, 'logfile' : cur_logfile, 'start' : False };
            index_cmds[faidx_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create reference faidx", 'outfile' : faidx_file, 'logfile' : cur_logfile, 'start' : False };

        if step == "index":
            cur_logfile = os.path.join(globs['iterlogdir'], "bwa-index-iter-" + globs['iter-str'] + ".log");
            index_files = [cur_ref + ".amb", cur_ref + ".ann", cur_ref + ".bwt", cur_ref + ".pac", cur_ref + ".sa"];

            index_cmd = globs['bwa-path'] + " index " + cur_ref;
            cmds[index_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create BWA reference index", 'outfile' : "", 'logfile' : cur_logfile, 'start' : False };
            index_cmds[index_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Create BWA reference index", 'outfile' : "", 'logfile' : cur_logfile, 'start' : False };

    index_procs = min(3, globs['num-procs']);
    pool = mp.Pool(processes=index_procs);
    for result in pool.starmap(PC.runCMD, ((index_cmd, globs, cmds, True) for index_cmd in index_cmds )):
        if result:
            pool.terminate();
            globs['exit-code'] = 1;
            PC.endProg(globs);
    pool.terminate();

    return cmds;

#############################################################################