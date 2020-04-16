# Mapping modules for Pseudo-it
#############################################################################
import os, lib.picore as PC
#############################################################################

def BWA(lib_item):
# Map a set of reads with BWA mem.
    cur_lib_type, cur_lib_files, cur_ref, globs, step_start_time = lib_item;

    cur_logfile = os.path.join(globs['logdir'], "bwa-mem-" + cur_lib_type + "-iter-" + globs['iter-str'] + ".log");
    bamfile = os.path.join(globs['iterdir'], cur_lib_type + "-iter-" + globs['iter-str'] + ".bam.gz");

    run_flag = PC.runCheck([bamfile], cur_logfile, globs);

    if run_flag:
        bwa_cmd = globs['bwa-path'] + " mem -t " + str(globs['bwa-t']) + " " + cur_ref + " " + cur_lib_files + " | " + globs['samtools-path'] + " view -bh - > " + bamfile;
        exit_flag = PC.runCMD(bwa_cmd, "BWA mem", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# " + cur_lib_type + " BAM file already exists", globs['pad'], sep=".") + bamfile);
        exit_flag = False;

    return bamfile, exit_flag;

#############################################################################

def mergeBam(bamfiles, globs):
# Merge BAM files from different library types.
    cur_logfile = os.path.join(globs['logdir'], "picard-merge-sam-iter-" + globs['iter-str'] + ".log");
    merged_bamfile = os.path.join(globs['iterdir'], "merged-iter-" + globs['iter-str'] + ".bam.gz");

    run_flag = PC.runCheck([merged_bamfile], cur_logfile, globs);

    if run_flag:
        merge_cmd = globs['picard-path'] + " MergeSamFiles ";
        for bam in bamfiles:
            merge_cmd += "I=" + bam + " ";
        merge_cmd += "USE_THREADING=TRUE VALIDATION_STRINGENCY=LENIENT O=" + merged_bamfile;
        exit_flag = PC.runCMD(merge_cmd, "BWA merge", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Merged BAM file already exists", globs['pad'], sep=".") + merged_bamfile + "\n");
        exit_flag = False;

    return merged_bamfile, exit_flag;

#############################################################################