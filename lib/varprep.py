# Pre-processing steps for variant calling
#############################################################################
import os, lib.picore as PC
#############################################################################

def addRG(bamfile, globs):
# Run Picard's AddOrReplaceReadGroups on a merged BAM file.
    cur_logfile = os.path.join(globs['logdir'], "picard-add-rg-iter-" + globs['iter-str'] + ".log");
    rg_bamfile = os.path.join(globs['iterdir'], "merged-rg-iter-" + globs['iter-str'] + ".bam.gz");
    rg_lib = "rg-iter-" + globs['iter-str'];

    run_flag = PC.runCheck([rg_bamfile], cur_logfile, globs);

    if run_flag:
        rg_cmd = globs['picard-path'] + " AddOrReplaceReadGroups I=" + bamfile + " O=" + rg_bamfile + " SO=coordinate LB=" + rg_lib + " PL=illumina PU=misc SM=" + rg_lib + " VALIDATION_STRINGENCY=LENIENT";
        if globs['tmpdir'] != "System default.":
            rg_cmd += " TMP_DIR=\"" + globs['tmpdir'] + "\"";
        exit_flag = PC.runCMD(rg_cmd, "Picard AddOrReplaceReadGroups", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# BAM with read groups already exists", globs['pad'], sep=".") + rg_bamfile + "\n");
        exit_flag = False;
    return rg_bamfile, exit_flag;

#############################################################################

def markDups(bamfile, globs):
# Mark duplicates of a BAM file
    cur_logfile = os.path.join(globs['logdir'], "picard-mkdup-iter-" + globs['iter-str'] + ".log");
    mkdup_bamfile = os.path.join(globs['iterdir'], "merged-rg-mkdup-iter-" + globs['iter-str'] + ".bam.gz");
    dupmet_file = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-dupmets.txt");

    run_flag = PC.runCheck([mkdup_bamfile], cur_logfile, globs);

    if run_flag:
        mkdup_cmd = globs['picard-path'] + " MarkDuplicates I=" + bamfile + " O=" + mkdup_bamfile + " VALIDATION_STRINGENCY=LENIENT M=" + dupmet_file;
        if globs['tmpdir'] != "System default.":
            mkdup_cmd += " TMP_DIR=\"" + globs['tmpdir'] + "\"";
        exit_flag = PC.runCMD(mkdup_cmd, "Picard MarkDuplicates", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# BAM with marked duplicates already exists", globs['pad'], sep=".") + mkdup_bamfile + "\n");
        exit_flag = False;
    return mkdup_bamfile, exit_flag;

#############################################################################

def indexBAM(bamfile, globs):
# Index a BAM file with samtools.
    cur_logfile = os.path.join(globs['logdir'], "samtools-index-iter-" + globs['iter-str'] + ".log");
    index_bamfile = bamfile + ".bai";

    run_flag = PC.runCheck([index_bamfile], cur_logfile, globs);

    if run_flag:
        index_cmd = globs['samtools-path'] + " index " + bamfile;
        exit_flag = PC.runCMD(index_cmd, "samtools index", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# BAM index already exists", globs['pad'], sep=".") + index_bamfile + "\n");
        exit_flag = False;
    return exit_flag;

#############################################################################
