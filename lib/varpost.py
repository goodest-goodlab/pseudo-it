# Variant post-processing for Pseudo-it
#############################################################################
import os, multiprocessing as mp, lib.picore as PC, lib.piref as piref
#############################################################################

def selectSNPs(globs, cmds, vcf_file):
# Run the command to select only SNPs from a VCF file.

    gatk_cmd = globs['gatk-path'] + " SelectVariants -V " + vcf_file + " -O " + globs['iter-final-vcf'] + " -select-type SNP -xl-select-type INDEL -xl-select-type MIXED -xl-select-type SYMBOLIC";
    cmd_num = PC.getCMDNum(globs, len(cmds));
    cmds[gatk_cmd] = { 'cmd-num' : cmd_num, 'desc' : "Select SNPs", 'outfile' : vcf_file,  'logfile' : globs['iter-final-vcf-log'], 'start' : False };
    exit_flag = PC.runCMD(gatk_cmd, globs, cmds, True);
    PC.exitCheck(exit_flag, globs);
    # End the program if an error is encountered

    return cmds;

#############################################################################

def varFilter(globs, cmds, cur_ref):
# Run the command to filter variants from a VCF file based on input filters. Default: "MQ < 30.0 || DP < 5 || DP > 60"
    
    bcftools_cmds = {};
    for scaff in globs['scaffolds']:
        # if not globs['last-iter'] or (globs['last-iter'] and not globs['indels']):
        #     cur_logfile = os.path.join(globs['itervcflogdir'], "bcftools-filter-" + scaff + "-iter-" + globs['iter-str'] + "-snps.log");
        #     vcf_file = os.path.join(globs['itervcfscaffdir'], scaff + "-iter-" + globs['iter-str'] + "-snps.vcf.gz");
        #     filter_file = os.path.join(globs['itervcfscaffdir'], scaff + "-iter-" + globs['iter-str'] + "-snps-filter.vcf.gz");
        # else:
        cur_logfile = os.path.join(globs['itervcflogdir'], "bcftools-filter-" + scaff + "-iter-" + globs['iter-str'] + ".log");
        vcf_file = os.path.join(globs['itervcfscaffdir'], scaff + "-iter-" + globs['iter-str'] + ".vcf.gz");
        filter_file = os.path.join(globs['itervcfscaffdir'], scaff + "-iter-" + globs['iter-str'] + "-filter.vcf.gz");

        bcftools_cmd = globs['bcftools-path'] + " filter -m+ -e " + globs['filter'] + "  -e 'ALT=\"*\"' -s pseudoit --IndelGap 5 -Oz -o " + filter_file + " " + vcf_file;

        cmd_num = PC.getCMDNum(globs, len(cmds));
        cmds[bcftools_cmd] = { 'cmd-num' : cmd_num, 'desc' : "Filter VCF " + scaff, 'outfile' : filter_file,  'logfile' : cur_logfile, 'start' : False, "vcffile" : vcf_file };
        bcftools_cmds[bcftools_cmd] = { 'cmd-num' : cmd_num, 'desc' : "Filter VCF " + scaff, 'outfile' : filter_file,  'logfile' : cur_logfile, 'start' : False, "vcffile" : vcf_file };

    if globs['dryrun']:
        cmd_num = PC.getCMDNum(globs, len(cmds));
        bcftools_skeleton_cmd = globs['bcftools-path'] + " filter -m+ -e " + globs['filter'] + " -s pseudoit --IndelGap 5 -Oz -o <filtered vcf> <input vcf>";
        cmds[bcftools_skeleton_cmd] = { 'cmd-num' : cmd_num, 'desc' : str(globs['num-procs']) + " bcftools filter procs in parallel", 'outfile' : "",  'logfile' : "", 'start' : False };
        PC.report_step(globs, cmds, bcftools_skeleton_cmd, "DRYRUN", bcftools_skeleton_cmd);

    else:
        pool = mp.Pool(processes=globs['filter-procs']);
        for result in pool.starmap(PC.runCMD, ((bcftools_cmd, globs, cmds, True) for bcftools_cmd in bcftools_cmds )):
            if result:
                pool.terminate();
                globs['exit-code'] = 1;
                PC.endProg(globs);
        pool.terminate();    

    return cmds;

#############################################################################

def gatherVCFs(globs, cmds):
# Combine the region VCFs from haplotypeCallerMulti.    

    # vcf_file = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + "-filter.vcf.gz");
    # cur_logfile = os.path.join(globs['iterlogdir'], "gatk-gathervcfs-iter-" + globs['iter-str'] + ".log");
    # if globs['last-iter'] and globs['indels']:
    #     vcf_file = vcf_file.replace("-filter.vcf.gz", "-filter-final.vcf.gz")
    #     cur_logfile = cur_logfile.replace(".log", "-final.log");
    params_file = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + "-gathervcfs-params.txt");        

    # infile_ext = "-snps-filter.vcf.gz";
    # if globs['last-iter']:
    #     if globs['indels']:
    #         infile_ext = "-filter.vcf.gz";
    infile_ext = "-filter.vcf.gz";

    with open(params_file, "w") as paramsfile:
        for scaff in globs['scaffolds']:
            scaff_vcf = os.path.join(globs['itervcfscaffdir'], scaff + "-iter-" + globs['iter-str'] + infile_ext);
            paramsfile.write("-I " + scaff_vcf + "\n");
    gatk_cmd = globs['gatk-path'] + " GatherVcfs --arguments_file " + params_file + " -O " + globs['iter-gather-vcf'];
    cmds[gatk_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Gather VCFs", 'outfile' : globs['iter-gather-vcf'],  'logfile' : globs['iter-gather-vcf-log'], 'start' : False };

    exit_flag = PC.runCMD(gatk_cmd, globs, cmds, True);
    PC.exitCheck(exit_flag, globs);
    # End the program if an error is encountered

    return cmds;
#############################################################################

def indexVCF(globs, cmds, vcf_file):
    
    index_file = vcf_file + ".tbi";
    cur_logfile = cur_logfile = os.path.join(globs['iterlogdir'], "tabix-" + globs['iter-str'] + ".log");

    index_cmd = "tabix -fp vcf " + vcf_file;
    cmds[index_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Index VCF", 'outfile' : index_file,  'logfile' : cur_logfile, 'start' : False };
    exit_flag = PC.runCMD(index_cmd, globs, cmds, True);
    PC.exitCheck(exit_flag, globs);
    # End the program if an error is encountered    

    return cmds;

#############################################################################










