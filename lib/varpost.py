# Variant post-processing for Pseudo-it
#############################################################################
import os, lib.picore as PC, lib.piref as piref
#############################################################################

def selectSNPs(vcffile, cur_ref, globs):
# Run the command to select only SNPs from a VCF file.

    cur_logfile = os.path.join(globs['logdir'], "gatk-selectvariants-iter-" + globs['iter-str'] + ".log");
    snp_vcffile = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-snps.vcf.gz");

    run_flag = PC.runCheck([snp_vcffile], cur_logfile, globs);

    if run_flag:
        gatk_cmd = globs['gatk-path'] + " SelectVariants -R " + cur_ref + " -V " + vcffile + " -O " + snp_vcffile + " -select-type SNP";
        exit_flag = PC.runCMD(gatk_cmd, "GATK SelectVariants", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# VCF file already exists", globs['pad'], sep=".") + snp_vcffile + "\n");
        exit_flag = False;

    return snp_vcffile, exit_flag;

#############################################################################

def varFilter(filter_item):
# Run the command to filter variants from a VCF file based on input filters. Default: "MQ < 30.0 || DP < 5 || DP > 60"
    vcffile, vcfdir, logdir, cur_ref, globs = filter_item;
    indel_suffix = "";
    if globs['iteration'] == globs['num-iters'] and not globs['indels']:
        indel_suffix = "-noindels"
    
    if globs['iteration'] != globs['num-iters'] or globs['num-procs'] == 1:
        cur_logfile = os.path.join(globs['logdir'], "gatk-varfilter-snps-iter-" + globs['iter-str'] + indel_suffix + ".log");
        filter_vcffile = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-filtered" + indel_suffix + ".vcf.gz");
    else:
        scaff = os.path.basename(vcffile);
        print("AHHHHHHHHHHH:" + scaff);
        scaff = scaff[:scaff.index("-iter")];
        cur_logfile = os.path.join(logdir, "gatk-varfilter-" + scaff + "-iter-" + globs['iter-str'] + ".log");
        filter_vcffile = os.path.join(vcfdir, scaff + "-iter-" + globs['iter-str'] + "-filtered.vcf.gz");
        
    run_flag = PC.runCheck([filter_vcffile], cur_logfile, globs);

    if run_flag:
        #gatk_cmd = globs['gatk-path'] + " VariantFiltration -R " + cur_ref + " -V " + vcffile + " -filter " + globs['filter'] + " -filter-name 'pseudoit-filter' -O " + filter_vcffile;
        bcftools_cmd = globs['bcftools-path'] + " filter -m+ -e " + globs['filter'] + " -s pseudoit --IndelGap 5 -Oz -o " + filter_vcffile + " " + vcffile;
        
        # PC.printWrite(globs['logfilename'], globs['log-v'], gatk_cmd);
        # with open(cur_logfile, "w") as lf:
        #     cmd_result = subprocess.run(gatk_cmd, shell=True, stdout=lf, stderr=lf);
        # #p = subprocess.Popen(my_cmd, shell=True)
        # #os.waitpid(p.pid, 0)
        # exit_flag = False;

        exit_flag = PC.runCMD(bcftools_cmd, "bcftools filter", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# VCF file already exists", globs['pad'], sep=".") + filter_vcffile + "\n");
        exit_flag = False;

    return filter_vcffile, exit_flag;

#############################################################################

def genConsensus(vcffile, cur_ref, globs):
# Run the command to generate a consensus FASTA file from the reference and the variants.
    cur_logfile = os.path.join(globs['logdir'], "gatk-genconsensus-iter-" + globs['iter-str'] + ".log");
    consensus_file = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + ".fa");

    run_flag = PC.runCheck([consensus_file], cur_logfile, globs);

    if run_flag:
        gatk_cmd = globs['gatk-path'] + " FastaAlternateReferenceMaker -R " + cur_ref + " -V " + vcffile + " -O " + consensus_file;
        exit_flag = PC.runCMD(gatk_cmd, "GATK FastaAlternateReferenceMaker", cur_logfile, True, globs);
        #bcftools_cmd = globs['bcftools-path'] + " consensus -f " + cur_ref + " -o " + consensus_file + " -e \"FILTER='pseudoit-filter'\" " + vcffile
        #exit_flag = PC.runCMD(bcftools_cmd, "bcftools consensus", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Consensus FASTA already exists", globs['pad'], sep=".") + consensus_file + "\n");
        exit_flag = False;

    return consensus_file, exit_flag;

#############################################################################

def fixHeaders(fafile, globs):
# Fix the headers from the consensus FASTA file.
    fa_fixed = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-final.fa");
    sed_file = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-fixheaders.sed");
    exit_flag = False;

    if not os.path.isfile(fa_fixed) or globs['overwrite']:
        # cur_scaffs, null = piref.getScaffs(fafile, globs, None, False);
        # scaff_rep = {};
        # for scaff in cur_scaffs:
        #     scaff_orig = scaff[scaff.index(" ")+1:scaff.rindex(":")];
        #     scaff_rep[scaff] = scaff_orig;


        # with open(sed_file, "w") as sf:
        #     for scaff in scaff_rep:
        #         sf.write("s/" + scaff + "/" + scaff_rep[scaff] + "/g;\n");

        # cmd = "sed -f " + sed_file + " < " + fafile + " > " + fa_fixed;

        cmd = "sed 's/>[0-9]* />/g; s/:1-[0-9]*//g' " + fafile + " > " + fa_fixed;

        if globs['dryrun']:
            PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : " + cmd);
        else:
            PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : " + cmd);
            os.system(cmd);

    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Fixed FASTA already exists", globs['pad'], sep=".") + fa_fixed + "\n");

    return fa_fixed, exit_flag;

#############################################################################

def getMask(vcffile, globs):
# Get the sites to be masked into a bed file.
    maskfile = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-masksites.bed");
    exit_flag = False;

    if not os.path.isfile(maskfile) or globs['overwrite']:
        #cmd = "zgrep \"\./\.\" " + vcffile + " > " + maskfile;
        cmd = "zgrep \"\./\.\" " + vcffile + " | awk '{{OFS=\"\t\"; if ($0 !~ /\#/); print $1, $2-1, $2}}' | bedtools merge -i - > " + maskfile;
        if globs['dryrun']:
            PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : " + cmd);
        else:
            PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : " + cmd);
            os.system(cmd);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Mask sites file already exists", globs['pad'], sep=".") + maskfile + "\n");

    return maskfile, exit_flag;

#############################################################################

def maskFa(maskbed, cur_ref, globs):
# Fix the headers from the consensus FASTA file.
    cur_logfile = os.path.join(globs['logdir'], "bedtools-maskfasta-" + globs['iter-str'] + ".log");
    mask_fa = os.path.join(globs['iterdir'], "iter-" + str(int(globs['iter-str']) - 1) + "-masked.fa");

    run_flag = PC.runCheck([mask_fa], cur_logfile, globs);

    if run_flag:
        mask_cmd = globs['bedtools-path'] + " maskfasta -fi " + cur_ref + " -bed " + maskbed + " -soft -fo " + mask_fa;
        exit_flag = PC.runCMD(mask_cmd, "tabix", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Masked FASTA file already exists", globs['pad'], sep=".") + mask_fa + "\n");
        exit_flag = False;

    return mask_fa, exit_flag;

#############################################################################

def rmNonRef(vcffile, globs):
# Deals with the weird <NON_REF> alleles in gatk's new GVCF format.
    vcffile_decom = vcffile.replace(".gz", "");

    cmd = "gunzip " + vcffile;
    if globs['dryrun']:
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : " + cmd);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : " + cmd);
        os.system(cmd);

    cmd = "sed -i 's/,<NON_REF>//g;s/<NON_REF>/./g;' " + vcffile_decom;
    if globs['dryrun']:
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : " + cmd);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : " + cmd);
        os.system(cmd);

    cmd = "bgzip " + vcffile_decom;
    if globs['dryrun']:
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : " + cmd);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : " + cmd);
        os.system(cmd);

#############################################################################

def finalConsensus(vcffile, cur_ref, globs):
# Run the command to generate a consensus FASTA file from the reference and the variants.
    if globs['indels']:
        suffix = "-indels";
    else:
        suffix = "-noindels"

    cur_logfile = os.path.join(globs['logdir'], "bcftools-consensus-iter-" + globs['iter-str'] + suffix + ".log");
    chain_file = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + suffix + "-final.chain");
    consensus_file = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + suffix + "-final.fa");

    if not globs['dryrun']:
        first_lower, linestr_orig, linestr_repl = getConsCase(cur_ref);
        # This first_lower stuff is a hack to deal with bcftools consensus using the case of the first base in the reference fasta to inject variants.
        # Possibly resolved: https://github.com/samtools/bcftools/issues/1150#issuecomment-582407490
        # Need to test and make sure it is in official release before I remove this hack.
    else:
        first_lower, linestr_orig, linestr_repl = True, "a", "A";

    if first_lower:
        cmd = "sed -i '2 s/" + linestr_orig + "/" + linestr_repl + "/g' " + cur_ref;
        if globs['dryrun']:
            PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : sed to convert case of first ref char." );
        else:
            PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : sed to convert case of first ref char." );
            os.system(cmd);
    # Part of first_lower hack.

    run_flag = PC.runCheck([consensus_file], cur_logfile, globs);

    if run_flag:
        bcftools_cmd = globs['bcftools-path'] + " consensus -f " + cur_ref + " -o " + consensus_file
        if globs['indels']:
            bcftools_cmd += " -c " + chain_file
        bcftools_cmd += " -e \"FILTER='pseudoit'\" " + vcffile;
        exit_flag = PC.runCMD(bcftools_cmd, "bcftools consensus", cur_logfile, True, globs);
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Consensus FASTA already exists", globs['pad'], sep=".") + consensus_file + "\n");
        exit_flag = False;
    # Consensus command

    if first_lower:
        if not globs['dryrun']:
            first_lower, linestr_orig, linestr_repl = getConsCase(cur_ref);

        cmd = "sed -i '2 s/" + linestr_orig + "/" + linestr_repl + "/g' " + cur_ref;
        if globs['dryrun']:
            PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : sed to revert case of first consensus char.");
        else:
            PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : sed to revert case of first consensus char.");
            os.system(cmd);        

        if not globs['dryrun']:
            first_lower, linestr_orig, linestr_repl = getConsCase(consensus_file);

        cmd = "sed -i '2 s/" + linestr_orig + "/" + linestr_repl + "/g' " + consensus_file;
        if globs['dryrun']:
            PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : sed to revert case of first ref char.");
        else:
            PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : sed to revert case of first ref char.");
            os.system(cmd);
    # Part of first_lower hack.

    return consensus_file, exit_flag;

#############################################################################

def getConsCase(cur_ref):
# bcftools consensus uses the case of the first base in the reference file for injecting
# variants. This is a hack to make sure it is always upper case.

    first_lower = False;
    lines = 0;
    for line in open(cur_ref):
        line = line.strip();
        lines += 1;
        if lines == 2:
            linestr_orig = line;
            linelist = list(linestr_orig);
            if linelist[0].islower():
                linelist[0] = linelist[0].upper();
                first_lower = True;
            elif linelist[0].isupper():
                linelist[0] = linelist[0].lower();
            linestr_repl = "".join(linelist);
            break;

    return first_lower, linestr_orig, linestr_repl;

#############################################################################