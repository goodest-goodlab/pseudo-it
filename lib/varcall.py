
# Variant calling for Pseudo-it
#############################################################################
import os, math, multiprocessing as mp, lib.picore as PC
#############################################################################

def haplotypeCaller(globs, cmds, cur_ref, dup_bamfile):
# Run HaplotypeCaller for each scaffold.

    gatk_cmds = {};
    for scaff in globs['scaffolds']:
        cur_logfile = os.path.join(globs['itervcflogdir'], "gatk-haplotypcaller-" + scaff + "-iter-" + globs['iter-str'] + ".log");
        if globs['last-iter']:
            vcffile = os.path.join(globs['itervcfscaffdir'], scaff + "-iter-" + globs['iter-str'] + ".gvcf.gz");
        else:
            vcffile = os.path.join(globs['itervcfscaffdir'], scaff + "-iter-" + globs['iter-str'] + ".vcf.gz");
        
        gatk_cmd = globs['gatk-path'] + " HaplotypeCaller -R " + cur_ref + " -I " + dup_bamfile + " -L \"" + scaff + "\" -stand-call-conf 30 --native-pair-hmm-threads " + str(globs['gatk-t']);
        if globs['last-iter']:
            gatk_cmd += " -ERC GVCF";
        # The final iteration outputs GVCFs to properly emit all sites
        gatk_cmd += " -O " + vcffile;

        cmd_num = PC.getCMDNum(globs, len(cmds));
        cmds[gatk_cmd] = { 'cmd-num' : cmd_num, 'desc' : "HaplotypeCaller " + scaff, 'outfile' : vcffile,  'logfile' : cur_logfile, 'start' : False };
        gatk_cmds[gatk_cmd] = { 'cmd-num' : cmd_num, 'desc' : "HaplotypeCaller " + scaff, 'outfile' : vcffile,  'logfile' : cur_logfile, 'start' : False };

    if globs['dryrun']:
        cmd_num = PC.getCMDNum(globs, len(cmds));
        gatk_skeleton_cmd = globs['gatk-path'] + " HaplotypeCaller -R <reference fasta> -I <BAM file> -L \"<scaffold>\" -stand-call-conf 30 --native-pair-hmm-threads " + str(globs['gatk-t']);
        if globs['last-iter']:
            gatk_skeleton_cmd += " -ERC GVCF";
        # The final iteration outputs GVCFs to properly emit all sites
        gatk_skeleton_cmd += " -O <vcf file>";
        cmds[gatk_skeleton_cmd] = { 'cmd-num' : cmd_num, 'desc' : str(globs['gatk-procs']) + " HaplotypeCaller procs in parallel", 'outfile' : "",  'logfile' : "", 'start' : False };
        PC.report_step(globs, cmds, gatk_skeleton_cmd, "DRYRUN", gatk_skeleton_cmd);

    else:
        pool = mp.Pool(processes=globs['gatk-procs']);
        for exit_flag in pool.starmap(PC.runCMD, ((gatk_cmd, globs, cmds, True) for gatk_cmd in gatk_cmds )):
            if exit_flag:
                pool.terminate();
                globs['exit-code'] = 1;
                PC.endProg(globs);
        pool.terminate();    

    return cmds;

#############################################################################

def genotypeGVCFs(globs, cmds, cur_ref):
# Genotype the GVCFs from the last iteration by scaffold.

    gatk_cmds = {};
    for scaff in globs['scaffolds']:
        cur_logfile = os.path.join(globs['itervcflogdir'], "gatk-genotypegvcfs- " + scaff + "-iter-" + globs['iter-str'] + ".log");
        gvcf_file = os.path.join(globs['itervcfscaffdir'], scaff + "-iter-" + globs['iter-str'] + ".gvcf.gz");
        vcf_file = os.path.join(globs['itervcfscaffdir'], scaff + "-iter-" + globs['iter-str'] + ".vcf.gz");

        gatk_cmd = globs['gatk-path'] + " GenotypeGVCFs -R " + cur_ref + " -V " + gvcf_file + " -O " + vcf_file + " --include-non-variant-sites";  

        cmd_num = PC.getCMDNum(globs, len(cmds));
        cmds[gatk_cmd] = { 'cmd-num' : cmd_num, 'desc' : "Genotype gVCF " + scaff, 'outfile' : vcf_file,  'logfile' : cur_logfile, 'start' : False };
        gatk_cmds[gatk_cmd] = { 'cmd-num' : cmd_num, 'desc' : "Genotype gVCF " + scaff, 'outfile' : vcf_file,  'logfile' : cur_logfile, 'start' : False };

    pool = mp.Pool(processes=globs['num-procs']);
    for exit_flag in pool.starmap(PC.runCMD, ((gatk_cmd, globs, cmds, True) for gatk_cmd in gatk_cmds )):
        if exit_flag:
            pool.terminate();
            globs['exit-code'] = 1;
            PC.endProg(globs);
    pool.terminate();    

    return cmds;


# #############################################################################

# def gatherVcfs(vcfdir, cur_ref, globs):
# # Combine the region VCFs from haplotypeCallerMulti.
#     cur_logfile = os.path.join(globs['logdir'], "gatk-gathervcfs-iter-" + globs['iter-str'] + ".log");
#     infile_ext = ".vcf.gz";
#     outfile_ext = ".vcf.gz";
#     if globs['iteration'] == globs['num-iters']:
#         infile_ext = "-filtered.vcf.gz";
#         outfile_ext = "-filtered-final.vcf.gz";
#         cur_logfile = cur_logfile.replace(".log", "-final.log");
#     vcf_file = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + outfile_ext);

#     run_flag = PC.runCheck([vcffile], cur_logfile, globs);

#     if run_flag:
#         params_file = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-gathervcfs-params.txt");
#         with open(params_file, "w") as paramsfile:
#             for scaff in globs['scaffolds']:
#                 scaff_vcf = os.path.join(vcfdir, scaff + "-iter-" + globs['iter-str'] + infile_ext);
#                 paramsfile.write("-I " + scaff_vcf + "\n");
#         gatk_cmd = globs['gatk-path'] + " GatherVcfs --arguments_file " + params_file + " -O " + vcffile;
#         exit_flag = PC.runCMD(gatk_cmd, "GATK GatherVcfs", cur_logfile, True, globs);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# VCF file already exists", globs['pad'], sep=".") + vcffile + "\n");
#         exit_flag = False;

#     return vcffile, exit_flag;

# #############################################################################

# def indexVCF(vcffile, globs, suffix=""):
# # Index the combined VCF from gatherVcfs.
#     if suffix != "":
#         suffix = "-" + suffix;
#     cur_logfile = os.path.join(globs['iterlogdir'], "vcf-index-iter-" + globs['iter-str'] + suffix + ".log");
#     if globs['iteration'] == globs['num-iters']:
#         cur_logfile = cur_logfile.replace(".log", "-final.log");
#     index_file = vcffile + ".tbi";

#     run_flag = PC.runCheck([index_file], cur_logfile, globs);

#     if run_flag:
#         index_cmd = "tabix -fp vcf " + vcffile;
#         exit_flag = PC.runCMD(index_cmd, "tabix", cur_logfile, True, globs);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# VCF index file already exists", globs['pad'], sep=".") + vcffile + "\n");
#         exit_flag = False;

#     return exit_flag;
    
# #############################################################################