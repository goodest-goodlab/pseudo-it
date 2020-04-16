
# Variant calling for Pseudo-it
#############################################################################
import os, lib.picore as PC
#############################################################################

def haplotypeCaller(bamfile, cur_ref, globs):
# Run HaplotypeCaller for a whole sequence.
    cur_logfile = os.path.join(globs['logdir'], "gatk-haplotypcaller-iter-" + globs['iter-str'] + ".log");
    
    if globs['iteration'] != globs['num-iters']:
        vcffile = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + ".vcf.gz");
    else:
        vcffile = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + ".gvcf.gz");

    run_flag = PC.runCheck([vcffile], cur_logfile, globs);

    if run_flag:
        gatk_cmd = globs['gatk-path'] + " HaplotypeCaller -R " + cur_ref + " -I " + bamfile + " -stand-call-conf 30 --native-pair-hmm-threads " + str(globs['gatk-t']);
        if globs['iteration'] == globs['num-iters']:
            gatk_cmd += " -ERC GVCF";
        # The final iteration outputs a GVCF to properly emit all sites.
        gatk_cmd += " -O " + vcffile;
        
        exit_flag = PC.runCMD(gatk_cmd, "GATK HaplotypeCaller", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# VCF file already exists:", globs['pad']) + vcffile + "\n");
        exit_flag = False;

    return vcffile, exit_flag;
#############################################################################

def haplotypeCallerMulti(vcf_item):
# Run HaplotypeCaller on a single region.
    scaff, bamfile, cur_ref, vcfdir, logdir, globs, step_start_time = vcf_item;
    #step_start_time = PC.report_stats(globs, "--> HaplotypeCaller : " + scaff, step_start=step_start_time);
    cur_logfile = os.path.join(logdir, "gatk-haplotypcaller-" + scaff + "-iter-" + globs['iter-str'] + ".log");
    
    if globs['iteration'] != globs['num-iters']:
        vcffile = os.path.join(vcfdir, scaff + "-iter-" + globs['iter-str'] + ".vcf.gz");
    else:
        vcffile = os.path.join(vcfdir, scaff + "-iter-" + globs['iter-str'] + ".gvcf.gz");

    run_flag = PC.runCheck([vcffile], cur_logfile, globs);

    if run_flag:
        gatk_cmd = globs['gatk-path'] + " HaplotypeCaller -R " + cur_ref + " -I " + bamfile + " -L \"" + scaff + "\" -stand-call-conf 30 --native-pair-hmm-threads " + str(globs['gatk-t']);
        if globs['iteration'] == globs['num-iters']:
            gatk_cmd += " -ERC GVCF";
        # The final iteration outputs GVCFs to properly emit all sites.
        gatk_cmd += " -O " + vcffile;
        exit_flag = PC.runCMD(gatk_cmd, "GATK HaplotypeCaller " + scaff, cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# VCF file already exists", globs['pad'], sep=".") + vcffile + "\n");
        exit_flag = False;

    return step_start_time, exit_flag;

#############################################################################

def genotypeGVCFs(gvcffile, cur_ref, globs):
# Run HaplotypeCaller for a whole sequence.
    cur_logfile = os.path.join(globs['logdir'], "gatk-genotypegvcfs-iter-" + globs['iter-str'] + ".log");
    vcffile = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + ".vcf.gz");

    run_flag = PC.runCheck([vcffile], cur_logfile, globs);

    if run_flag:
        gatk_cmd = globs['gatk-path'] + " GenotypeGVCFs -R " + cur_ref + " -V " + gvcffile + " -O " + vcffile + " --include-non-variant-sites";  
        exit_flag = PC.runCMD(gatk_cmd, "GATK GenotypeGVCFs", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# VCF file already exists:", globs['pad']) + vcffile + "\n");
        exit_flag = False;

    return vcffile, exit_flag;

#############################################################################

def genotypeGVCFsMulti(gvcf_item):
# Run HaplotypeCaller on a single region.
    scaff, cur_ref, vcfdir, logdir, globs, step_start_time = gvcf_item;
    #step_start_time = PC.report_stats(globs, "--> HaplotypeCaller : " + scaff, step_start=step_start_time);
    cur_logfile = os.path.join(logdir, "gatk-genotypegvfs-" + scaff + "-iter-" + globs['iter-str'] + ".log");
    gvcffile = os.path.join(vcfdir, scaff + "-iter-" + globs['iter-str'] + ".gvcf.gz");
    vcffile = os.path.join(vcfdir, scaff + "-iter-" + globs['iter-str'] + ".vcf.gz");

    run_flag = PC.runCheck([vcffile], cur_logfile, globs);

    if run_flag:
        gatk_cmd = globs['gatk-path'] + " GenotypeGVCFs -R " + cur_ref + " -V " + gvcffile + " -L \"" + scaff + "\" -O " + vcffile + " --include-non-variant-sites";
        exit_flag = PC.runCMD(gatk_cmd, "GATK HaplotypeCaller " + scaff, cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# VCF file already exists", globs['pad'], sep=".") + vcffile + "\n");
        exit_flag = False;

    return step_start_time, exit_flag;

#############################################################################

def gatherVcfs(vcfdir, cur_ref, globs):
# Combine the region VCFs from haplotypeCallerMulti.
    cur_logfile = os.path.join(globs['logdir'], "gatk-gathervcfs-iter-" + globs['iter-str'] + ".log");
    gather_ext = ".vcf.gz";
    if globs['iteration'] == globs['num-iters']:
        gather_ext = "-filtered.vcf.gz";
    vcffile = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + gather_ext);

    run_flag = PC.runCheck([vcffile], cur_logfile, globs);

    if run_flag:
        params_file = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-gathervcfs-params.txt");
        with open(params_file, "w") as paramsfile:
            for scaff in globs['scaffolds']:
                scaff_vcf = os.path.join(vcfdir, scaff + "-iter-" + globs['iter-str'] + gather_ext);
                paramsfile.write("-I " + scaff_vcf + "\n");
        gatk_cmd = globs['gatk-path'] + " GatherVcfs --arguments_file " + params_file + " -O " + vcffile;
        exit_flag = PC.runCMD(gatk_cmd, "GATK GatherVcfs", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# VCF file already exists", globs['pad'], sep=".") + vcffile + "\n");
        exit_flag = False;

    return vcffile, exit_flag;

#############################################################################

def indexVCF(vcffile, globs, suffix=""):
# Index the combined VCF from gatherVcfs.
    if suffix != "":
        suffix = "-" + suffix;
    cur_logfile = os.path.join(globs['logdir'], "vcf-index-iter-" + globs['iter-str'] + suffix + ".log");
    index_file = vcffile + ".tbi";

    run_flag = PC.runCheck([index_file], cur_logfile, globs);

    if run_flag:
        index_cmd = "tabix -fp vcf " + vcffile;
        exit_flag = PC.runCMD(index_cmd, "tabix", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# VCF index file already exists", globs['pad'], sep=".") + vcffile + "\n");
        exit_flag = False;

    return exit_flag;
    
#############################################################################