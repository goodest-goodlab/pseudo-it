# Mapping modules for Pseudo-it
#############################################################################
import os, gzip, math, multiprocessing as mp, lib.picore as PC
#############################################################################

def BWA(globs, cmds, cur_ref):
# Map a set of reads with BWA mem.

    bwa_cmds, bamfiles = {}, [];
    for lib_type in globs['libs']:
        if lib_type == "pe":
            fq_file = globs['libs'][lib_type].split(" ")[0];
        else:
            fq_file = globs['libs'][lib_type];
        # Get the fastq file for the current library.

        try:
            with gzip.open(fq_file) as fo:
                fl = fo.readline().decode().strip()[1:].split(":");
        except:
            with open(fq_file) as fo:
                fl = fo.readline().strip()[1:].split(":");
        # Try to read the first line whether the fastq file is compressed or not.        

        if len(fl) < 3:
            read_group_id = "pi-iter-" + globs['iter-str'];
        else:
            read_group_id = ".".join([fl[0], fl[1], fl[2]]);

        read_group = "@RG\\tID:" + read_group_id;
        read_group += "\\tPL:ILLUMINA";
        read_group += "\\tPU:" + read_group_id;
        read_group += "\\tLB:" + globs['sample-name'] + "-pseudoit-iter-" + globs['iter-str'] + "-" + lib_type;
        read_group += "\\tSM:" + globs['sample-name'];
        # Extract read group info from first line.
        # Assumes all reads in fastq file are from same run.
        # Assumes Illumina reads.
        # Sample and Library are determined here.

        cur_logfile = os.path.join(globs['iterlogdir'], "bwa-mem-" + lib_type + "-iter-" + globs['iter-str'] + ".log");
        bamfile = os.path.join(globs['iterbamdir'], lib_type + "-iter-" + globs['iter-str'] + ".bam.gz");
        bamfiles.append(bamfile);

        bwa_cmd = globs['bwa-path'] + " mem -t " + str(globs['bwa-t']) + " -R '" + read_group + "' " + cur_ref + " " + globs['libs'][lib_type] + " | " + globs['samtools-path'] + " sort | " + globs['samtools-path'] + " view -bh - > " + bamfile;
        cmd_num = PC.getCMDNum(globs, len(cmds))
        cmds[bwa_cmd] = { 'cmd-num' : cmd_num, 'desc' : "BWA " + lib_type + " read mapping", 'outfile' : bamfile,  'logfile' : cur_logfile, 'start' : False };

        bwa_cmds[bwa_cmd] = { 'cmd-num' : cmd_num, 'desc' : "BWA " + lib_type + " read mapping", 'logfile' : cur_logfile, 'start' : False };
    # Prepare the BWA commands for each library

    pool = mp.Pool(processes=globs['map-procs']);
    for result in pool.starmap(PC.runCMD, ((bwa_cmd, globs, cmds, True) for bwa_cmd in bwa_cmds )):
        if result:
            pool.terminate();
            globs['exit-code'] = 1;
            PC.endProg(globs);
    pool.terminate();
    # Run the BWA commands across multiple processors, if specified
    # End the program if an error is encountered

    return bamfiles, cmds;

#############################################################################

# def addRG(globs, cmds, bamfiles):
# # Run Picard's AddOrReplaceReadGroups on a merged BAM file.
# # This is now done during read mapping with bwa's -R option.

#     bwa_cmds, rg_bamfiles = {}, [];
#     for lib_type in bamfiles:
#         bamfile = bamfiles[lib_type];
#         cur_logfile = os.path.join(globs['iterlogdir'], "picard-add-rg-" + lib_type + "-iter-" + globs['iter-str'] + ".log");
#         rg_bamfile = os.path.join(globs['iterbamdir'], lib_type + "-iter-" + globs['iter-str'] + "-rg.bam.gz");
#         bamfiles.append(bamfile);

#         rg_cmd = globs['picard-path'] + " AddOrReplaceReadGroups I=" + bamfile + " O=" + rg_bamfile + " SO=coordinate LB=" + lib_type + " PL=illumina PU=misc SM=" + rg_lib + " VALIDATION_STRINGENCY=LENIENT";
#         if globs['tmpdir'] != "System default.":
#             rg_cmd += " TMP_DIR=\"" + globs['tmpdir'] + "\"";



#         bwa_cmd = globs['bwa-path'] + " mem -t " + str(globs['bwa-t']) + " " + cur_ref + " " + globs['libs'][lib_type] + " | " + globs['samtools-path'] + " view -bh - > " + bamfile;
#         cmd_num = PC.getCMDNum(globs, len(cmds))
#         cmds[bwa_cmd] = { 'cmd-num' : cmd_num, 'desc' : "BWA " + lib_type + " read mapping", 'outfile' : bamfile,  'logfile' : cur_logfile, 'start' : False };

#         bwa_cmds[bwa_cmd] = { 'cmd-num' : cmd_num, 'desc' : "BWA " + lib_type + " read mapping", 'logfile' : cur_logfile, 'start' : False };
#     # Prepare the BWA commands for each library

#     cur_logfile = os.path.join(globs['iterlogdir'], "picard-add-rg-iter-" + globs['iter-str'] + ".log");
#     rg_bamfile = os.path.join(globs['iterbamdir'], "merged-rg-iter-" + globs['iter-str'] + ".bam.gz");
#     rg_lib = "rg-iter-" + globs['iter-str'];

#     rg_cmd = globs['picard-path'] + " AddOrReplaceReadGroups I=" + merged_bamfile + " O=" + rg_bamfile + " SO=coordinate LB=" + rg_lib + " PL=illumina PU=misc SM=" + rg_lib + " VALIDATION_STRINGENCY=LENIENT";
#     if globs['tmpdir'] != "System default.":
#         rg_cmd += " TMP_DIR=\"" + globs['tmpdir'] + "\"";

#     cmds[rg_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Add read groups", 'outfile' : rg_bamfile, 'logfile' : cur_logfile, 'start' : False };
#     exit_flag = PC.runCMD(rg_cmd, globs, cmds, True);
#     PC.exitCheck(exit_flag, globs);

#     return rg_bamfile, cmds;

#############################################################################

def mergeBam(globs, cmds, bamfiles):
# Merge BAM files from different library types.

    cur_logfile = os.path.join(globs['iterlogdir'], "picard-merge-bam-iter-" + globs['iter-str'] + ".log");
    merged_bamfile = os.path.join(globs['iterbamdir'], "merged-iter-" + globs['iter-str'] + ".bam.gz");

    if len(bamfiles) > 1:
        merge_cmd = globs['picard-path'] + " MergeSamFiles ";
        for bamfile in bamfiles:
            merge_cmd += "I=" + bamfile + " ";
        if globs['tmpdir'] != "System default.":
            merge_cmd += "TMP_DIR=\"" + globs['tmpdir'] + "\" ";
        merge_cmd += "USE_THREADING=TRUE VALIDATION_STRINGENCY=LENIENT O=" + merged_bamfile;

        cmds[merge_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Merge BAM files", 'outfile' : merged_bamfile, 'logfile' : cur_logfile, 'start' : False };
        exit_flag = PC.runCMD(merge_cmd, globs, cmds, True);
        PC.exitCheck(exit_flag, globs);
        # End the program if an error is encountered

    else:
        merge_cmd = "mv " + bamfiles[0] + " " + merged_bamfile;
        cmds[merge_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Rename BAM file", 'outfile' : merged_bamfile, 'logfile' : "", 'start' : False };

        if globs['dryrun']:
            PC.report_step(globs, cmds, merge_cmd, "DRYRUN");
        else:
            PC.report_step(globs, cmds, merge_cmd, "EXECUTING");
            os.system(merge_cmd);
            if os.path.isfile(merged_bamfile):
                PC.report_step(globs, cmds, merge_cmd, "SUCCESS");
            else:
                PC.report_step(globs, cmds, merge_cmd, "ERROR");
                PC.errorOut("PIMAP1", "Error renaming BAM file.", globs);

    return merged_bamfile, cmds;

#############################################################################

def markDups(globs, cmds, rg_bamfile):
# Mark duplicates of a BAM file

    dupmet_file = os.path.join(globs['iterbamdir'], "iter-" + globs['iter-str'] + "-dupmets.txt");

    mkdup_cmd = globs['picard-path'] + " MarkDuplicates I=" + rg_bamfile + " O=" + globs['iter-final-bam'] + " VALIDATION_STRINGENCY=LENIENT M=" + dupmet_file + " CREATE_INDEX=true";
    if globs['tmpdir'] != "System default.":
        mkdup_cmd += " TMP_DIR=\"" + globs['tmpdir'] + "\"";

    cmds[mkdup_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Mark duplicates", 'outfile' : globs['iter-final-bam'], 'logfile' : globs['iter-final-bam-log'], 'start' : False };
    exit_flag = PC.runCMD(mkdup_cmd, globs, cmds, True);
    PC.exitCheck(exit_flag, globs);

    return cmds;

#############################################################################

# def indexBAM(globs, cmds):
# # Index a BAM file with samtools.
# # Now done durign MarkDuplicates with CREATE_INDEX=true

#     cur_logfile = os.path.join(globs['iterlogdir'], "samtools-index-iter-" + globs['iter-str'] + ".log");
#     index_bamfile = globs['iter-final-bam'] + ".bai";

#     index_cmd = globs['samtools-path'] + " index " + globs['iter-final-bam'];
#     cmds[index_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Index BAM file", 'outfile' : index_bamfile, 'logfile' : cur_logfile, 'start' : False };
#     exit_flag = PC.runCMD(index_cmd, globs, cmds, True);
#     PC.exitCheck(exit_flag, globs);

#     return cmds;

#############################################################################