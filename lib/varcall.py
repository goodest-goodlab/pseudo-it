
# Variant calling for Pseudo-it
#############################################################################
import os, math, multiprocessing as mp, lib.picore as PC
#############################################################################

def haplotypeCaller(globs, cmds, cur_ref, dup_bamfile):
# Run HaplotypeCaller for each specified interval.

    gatk_cmds = {};
    # A dictionary containing the gatk commands to be generated.

    for region in globs['intervals']:
        if globs['bed-mode'] == "file":
            region_str = "";
        else:
            region_str = region + "-";
        # For file names, if the region needs to be appended add a hyphen.

        cur_logfile = os.path.join(globs['itervcflogdir'], "gatk-haplotypcaller-" + region_str + "iter-" + globs['iter-str'] + ".log");
        if globs['last-iter'] and globs['mask'] != "none":
            vcffile = os.path.join(globs['itervcfscaffdir'], region_str + "iter-" + globs['iter-str'] + ".gvcf.gz");
        else:
            vcffile = os.path.join(globs['itervcfscaffdir'], region_str + "iter-" + globs['iter-str'] + ".vcf.gz");
        
        gatk_cmd = globs['gatk-path'] + " HaplotypeCaller -R " + cur_ref + " -I " + dup_bamfile + " -L \"" + region + "\" -stand-call-conf 30 --native-pair-hmm-threads " + str(globs['gatk-t']);
        if globs['last-iter'] and globs['mask'] != "none":
            gatk_cmd += " -ERC GVCF";
        # The final iteration outputs GVCFs to properly emit all sites if masking is to be done.
        gatk_cmd += " -O " + vcffile;

        cmd_num = PC.getCMDNum(globs, len(cmds));
        cmds[gatk_cmd] = { 'cmd-num' : cmd_num, 'desc' : "HaplotypeCaller " + region_str.replace("-",""), 'outfile' : vcffile,  'logfile' : cur_logfile, 'start' : False };
        gatk_cmds[gatk_cmd] = { 'cmd-num' : cmd_num, 'desc' : "HaplotypeCaller " + region_str.replace("-",""), 'outfile' : vcffile,  'logfile' : cur_logfile, 'start' : False };

    if globs['dryrun']:
        cmd_num = PC.getCMDNum(globs, len(cmds));
        gatk_skeleton_cmd = globs['gatk-path'] + " HaplotypeCaller -R <reference fasta> -I <BAM file> -L \"<region>\" -stand-call-conf 30 --native-pair-hmm-threads " + str(globs['gatk-t']);
        if globs['last-iter']:
            gatk_skeleton_cmd += " -ERC GVCF";
        # The final iteration outputs GVCFs to properly emit all sites if masking is to be done.
        gatk_skeleton_cmd += " -O <vcf file>";
        cmds[gatk_skeleton_cmd] = { 'cmd-num' : cmd_num, 'desc' : str(globs['gatk-procs']) + " HaplotypeCaller procs in parallel", 'outfile' : "",  'logfile' : "", 'start' : False };
        PC.report_step(globs, cmds, gatk_skeleton_cmd, "DRYRUN", gatk_skeleton_cmd);

    else:
        pool = mp.Pool(processes=globs['gatk-procs']);
        #pool_key = 'var-pool-' + str(globs['iteration']);
        #with pools[pool_key] as pool:
        for exit_flag in pool.imap_unordered(PC.runCMD, ((gatk_cmd, globs, cmds, True) for gatk_cmd in gatk_cmds )):
                if exit_flag:
                    pool.terminate();
                    globs['exit-code'] = 1;
                    PC.endProg(globs);
        pool.terminate();
    return cmds;

#############################################################################

def genotypeGVCFs(globs, cmds, cur_ref):
# Genotype the GVCFs from the last iteration by specified interval if masking is to be done.

    gatk_cmds = {};
    # A dictionary containing the gatk commands to be generated.

    for region in globs['intervals']:
        if globs['bed-mode'] == "file":
            region_str = "";
        else:
            region_str = region + "-";
        # For file names, if the region needs to be appended add a hyphen.

        cur_logfile = os.path.join(globs['itervcflogdir'], "gatk-genotypegvcfs- " + region_str + "iter-" + globs['iter-str'] + ".log");
        gvcf_file = os.path.join(globs['itervcfscaffdir'], region_str + "iter-" + globs['iter-str'] + ".gvcf.gz");
        vcf_file = os.path.join(globs['itervcfscaffdir'], region_str + "iter-" + globs['iter-str'] + ".vcf.gz");

        gatk_cmd = globs['gatk-path'] + " GenotypeGVCFs -R " + cur_ref + " -V " + gvcf_file + " -O " + vcf_file + " --include-non-variant-sites";  

        cmd_num = PC.getCMDNum(globs, len(cmds));
        cmds[gatk_cmd] = { 'cmd-num' : cmd_num, 'desc' : "Genotype gVCF " + region_str.replace("-",""), 'outfile' : vcf_file,  'logfile' : cur_logfile, 'start' : False };
        gatk_cmds[gatk_cmd] = { 'cmd-num' : cmd_num, 'desc' : "Genotype gVCF " + region_str.replace("-",""), 'outfile' : vcf_file,  'logfile' : cur_logfile, 'start' : False };


    if globs['dryrun']:
        cmd_num = PC.getCMDNum(globs, len(cmds));
        gatk_skeleton_cmd = globs['gatk-path'] + " GenotypeGVCFs -R <reference fasta> -V <gvcf file> -O <vcf file> --include-non-variant-sites";
        cmds[gatk_skeleton_cmd] = { 'cmd-num' : cmd_num, 'desc' : str(globs['gvcf-procs']) + " GenotypeGVCFs procs in parallel", 'outfile' : "",  'logfile' : "", 'start' : False };
        PC.report_step(globs, cmds, gatk_skeleton_cmd, "DRYRUN", gatk_skeleton_cmd);
    else:
        pool = mp.Pool(processes=globs['gvcf-procs']);
        for exit_flag in pool.imap_unordered(PC.runCMD, ((gatk_cmd, globs, cmds, True) for gatk_cmd in gatk_cmds )):
            if exit_flag:
                pool.terminate();
                globs['exit-code'] = 1;
                PC.endProg(globs);
        pool.terminate();    

    return cmds;

#############################################################################