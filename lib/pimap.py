# Mapping modules for Pseudo-it
#############################################################################
import os, gzip, math, multiprocessing as mp, lib.picore as PC
#############################################################################

def getRG(globs, cmds):
# This function reads the first line of a FASTQ file to get info
# for the read groups and stores this in the global dictionary to
# add to BAM files as reads are mapped.

    cmd = "getRG()";
    cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Getting read group information from FASTQ", 'outfile' : "", 'logfile' : "", 'start' : False };

    if not globs['dryrun']:
        PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
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

            globs['rg']['ID'] = read_group_id;
            globs['rg']['PL'] = "ILLUMINA";
            globs['rg']['PU'] = read_group_id;
            globs['rg']['LB'] = globs['sample-name'] + "-pseudoit-iter-" + globs['iter-str']# + "-" + lib_type;;
            globs['rg']['SM'] = globs['sample-name'];
            # Extract read group info from first line.
            # Assumes all reads in fastq file are from same run.
            # Assumes Illumina reads.
            # Sample and Library are determined here.

            break;
            # Only need to read the first file since we assume all are from the same library.

        PC.report_step(globs, cmds, cmd, "SUCCESS", "Read group info read", "");
    else:
        PC.report_step(globs, cmds, cmd, "DRYRUN", cmd);

    return globs, cmds;      

#############################################################################

def BWA(globs, cmds, cur_ref):
# Map a set of reads with BWA mem.

    bwa_cmds, bamfiles = {}, [];
    for lib_type in globs['libs']:
    # Generate a BWA command for each input fastq type.

        cur_logfile = os.path.join(globs['iterlogdir'], "bwa-mem-" + lib_type + "-iter-" + globs['iter-str'] + ".log");
        bamfile = os.path.join(globs['iterbamdir'], lib_type + "-iter-" + globs['iter-str'] + ".bam.gz");
        bamfiles.append(bamfile);
        # Get the bam file and log file for the current fastq file.

        rg_fields = ["ID", "PL", "PU", "LB", "SM"];
        rg_str = ["@RG"] + [ field + ":" + globs['rg'][field] for field in rg_fields ];
        rg_str = "\\t".join(rg_str);
        # Gets the read group info from globs and parses it for BWA's -R option

        bwa_cmd = globs['mapper-path'] + " mem -t " + str(globs['map-t']) + " -M -R '" + rg_str + "' " + cur_ref + " " + globs['libs'][lib_type];
        bwa_cmd += " | " + globs['samtools-path'] + " sort";
        bwa_cmd += " | " + globs['samtools-path'] + " view -bh -";
        bwa_cmd += " > " + bamfile;
        # Generate the bwa mem command for the current fastq file, including passing output to samtools for sorting and
        # converting to .bam.
        
        cmd_num = PC.getCMDNum(globs, len(cmds))
        # Get the current command number for the log.
        
        cmds[bwa_cmd] = { 'cmd-num' : cmd_num, 'desc' : "BWA " + lib_type + " read mapping", 'outfile' : bamfile,  'logfile' : cur_logfile, 'start' : False };
        # Save the bwa mem command to the global cmds dict.

        bwa_cmds[bwa_cmd] = { 'cmd-num' : cmd_num, 'desc' : "BWA " + lib_type + " read mapping", 'logfile' : cur_logfile, 'start' : False };
        # Save the bwa mem command to the bwa_cmds dict.
    # Prepare the BWA commands for each library

    pool = mp.Pool(processes=globs['map-procs']);
    for result in pool.imap_unordered(PC.runCMD, ((bwa_cmd, globs, cmds, True) for bwa_cmd in bwa_cmds )):
        if result:
            pool.terminate();
            globs['exit-code'] = 1;
            PC.endProg(globs);
    pool.terminate();
    # Run the BWA commands across multiple processors, if specified
    # End the program if an error is encountered

    return bamfiles, cmds;

#############################################################################

def hisat2(globs, cmds, cur_ref):
# Map a set of reads with BWA mem.

    hisat2_cmds, bamfiles = {}, [];
    for lib_type in globs['libs']:
    # Generate a hisat2 command for each input fastq type.

        cur_logfile = os.path.join(globs['iterlogdir'], "hisat2-" + lib_type + "-iter-" + globs['iter-str'] + ".log");
        bamfile = os.path.join(globs['iterbamdir'], lib_type + "-iter-" + globs['iter-str'] + ".bam.gz");
        bamfiles.append(bamfile);
        # Get the bam file and log file for the current fastq file.

        rg_fields = ["ID", "PL", "PU", "LB", "SM"];
        # The read group fields to add to the output bam.

        hisat2_cmd = globs['mapper-path'];
        for field in rg_fields:
            hisat2_cmd += " --rg " + field + ":" + globs['rg'][field];
        hisat2_cmd += " -p " + str(globs['map-t']);
        hisat2_cmd += " -x " + cur_ref;
        if lib_type == 'pe':
            hisat2_cmd += " -1 " + globs['libs'][lib_type].split(" ")[0];
            hisat2_cmd += " -2 " + globs['libs'][lib_type].split(" ")[1];
        else:
            hisat2_cmd += " -U " + globs['libs'][lib_type];
        hisat2_cmd += " | " + globs['samtools-path'] + " sort";
        hisat2_cmd += " | " + globs['samtools-path'] + " view -bh -";
        hisat2_cmd += " > " + bamfile;
        # Generate the hisat2 command, including adding read group info with --rg, and passing output to samtools for sorting
        # converting to .bam.

        cmd_num = PC.getCMDNum(globs, len(cmds));
        # Get the current command number for the log.

        cmds[hisat2_cmd] = { 'cmd-num' : cmd_num, 'desc' : "BWA " + lib_type + " read mapping", 'outfile' : bamfile,  'logfile' : cur_logfile, 'start' : False };
        # Save the hisat command to the global cmds dict.

        hisat2_cmds[hisat2_cmd] = { 'cmd-num' : cmd_num, 'desc' : "BWA " + lib_type + " read mapping", 'logfile' : cur_logfile, 'start' : False };
        # Save the hisat2 command to the bwa_cmds dict.
    # Prepare the hisat2 commands for each fastq type

    pool = mp.Pool(processes=globs['map-procs']);
    for result in pool.imap_unordered(PC.runCMD, ((hisat2_cmd, globs, cmds, True) for hisat2_cmd in hisat2_cmds )):
        if result:
            pool.terminate();
            globs['exit-code'] = 1;
            PC.endProg(globs);
    pool.terminate();
    # Run the hisat2 commands across multiple processors, if specified
    # End the program if an error is encountered

    return bamfiles, cmds;

#############################################################################

def mergeBam(globs, cmds, bamfiles):
# Merge BAM files from different library types.

    cur_logfile = os.path.join(globs['iterlogdir'], "picard-merge-bam-iter-" + globs['iter-str'] + ".log");
    merged_bamfile = os.path.join(globs['iterbamdir'], "merged-iter-" + globs['iter-str'] + ".bam.gz");
    # Get the log file and merged bam file name to output to.

    if len(bamfiles) > 1:
    # We only need to run picard if there are multiple bam files from mapping

        merge_cmd = globs['picard-path'] + " MergeSamFiles ";
        for bamfile in bamfiles:
            merge_cmd += "I=" + bamfile + " ";
        if globs['tmpdir'] != "System default.":
            merge_cmd += "TMP_DIR=\"" + globs['tmpdir'] + "\" ";
        if not globs['mkdups']:
            merge_cmd += "CREATE_INDEX=true ";
        merge_cmd += "USE_THREADING=TRUE VALIDATION_STRINGENCY=LENIENT O=" + merged_bamfile;
        # Generate the MergeSamFiles command.

        cmds[merge_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Merge BAM files", 'outfile' : merged_bamfile, 'logfile' : cur_logfile, 'start' : False };
        # Add the MergeSamFiles command to the global cmds dict.

        exit_flag = PC.runCMD((merge_cmd, globs, cmds, True));
        PC.exitCheck(exit_flag, globs);
        # Run the command and check for errors.

    else:
    # If there was only one bam file from mapping we don't need to merge, just move it to the expected location.

        merge_cmd = "mv " + bamfiles[0] + " " + merged_bamfile;
        # Generate the mv command.

        cmds[merge_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Rename BAM file", 'outfile' : merged_bamfile, 'logfile' : "", 'start' : False };
        # Add the mv command to the global commands dict.

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
        # Run the command and check for errors.

    return merged_bamfile, cmds;

#############################################################################

def markDups(globs, cmds, rg_bamfile):
# Mark duplicates of a BAM file

    dupmet_file = os.path.join(globs['iterbamdir'], "iter-" + globs['iter-str'] + "-dupmets.txt");
    # Get the duplicate metrics file name required to output by picard.

    mkdup_cmd = globs['picard-path'] + " MarkDuplicates I=" + rg_bamfile + " O=" + globs['iter-final-bam'] + " VALIDATION_STRINGENCY=LENIENT M=" + dupmet_file + " CREATE_INDEX=true";
    if globs['tmpdir'] != "System default.":
        mkdup_cmd += " TMP_DIR=\"" + globs['tmpdir'] + "\"";
    # Generate the MarkDuplicates command.

    cmds[mkdup_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Mark duplicates", 'outfile' : globs['iter-final-bam'], 'logfile' : globs['iter-final-bam-log'], 'start' : False };
    # Add the MarkDuplicates command to the global cmds dict.

    exit_flag = PC.runCMD((mkdup_cmd, globs, cmds, True));
    PC.exitCheck(exit_flag, globs);
    # Run the MarkDuplicates command and check for errors.

    return cmds;

#############################################################################
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

# def hisat2Index(globs, cmds, cur_ref):
# # Index a FASTA file with hisat2.
#     cur_logfile = os.path.join(globs['iterlogdir'], "hisat2-build-iter-" + globs['iter-str'] + ".log");
#     # Get the name of the logfile for the index command.

#     hisat2_build_cmd = globs['mapper-path'] + "-build " + cur_ref + " " + cur_ref;
#     # Generate the hisat-build command.

#     cmd_num = PC.getCMDNum(globs, len(cmds));
#     # Get the current command number for the log.

#     cmds[hisat2_build_cmd] = { 'cmd-num' : cmd_num, 'desc' : "hisat2-build index", 'outfile' : cur_ref,  'logfile' : cur_logfile, 'start' : False };
#     # Save the hisat-build command to the global cmds dict.

#     exit_flag = PC.runCMD(hisat2_build_cmd, globs, cmds, True);
#     PC.exitCheck(exit_flag, globs);
#     # Run the hisat-build command and check for errors.

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