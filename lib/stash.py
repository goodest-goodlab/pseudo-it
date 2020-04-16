#############################################################################
# Old or supplemental code that didn't end up being in the program

# WHEN I THOUGHT I NEEDED TO SPLIT OUT SNPS AND INDELS ON THE LAST ITERATION
    # if globs['iteration'] != globs['num-iters']:
    #     step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Select SNPs", step_start=step_start_time);
    #     vcf_snp_files = {};
    #     for result in pool.imap(varpost.selectDistributor, ((snp_type, vcf_file, cur_ref, globs) for snp_type in ['SNP', 'INDEL'])):
    #         snp_type, vcf_result, exit_flag = result;
    #         if exit_flag:
    #             pool.terminate();
    #             PC.endProg(globs);
    #         if vcf_result:
    #             vcf_snp_files[snp_type] = vcf_result;
    #     PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    #     # Select SNPs only, but also get indels on last iteration.

    #     step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Filter variants", step_start=step_start_time);
    #     vcf_filter_files = {};
    #     for result in pool.imap(varpost.varFilter, ((snp_type, vcf_snp_files[snp_type], cur_ref, globs) for snp_type in vcf_snp_files)):
    #         snp_type, vcf_result, exit_flag = result;
    #         if exit_flag:
    #             pool.terminate();
    #             PC.endProg(globs);
    #         vcf_filter_files[snp_type] = vcf_result;
    #     PC.exitCheck(exit_flag, globs);
    #     PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    #     # Filter VCF    




# def selectDistributor(select_item):
# # This function distributes the SelectVariants jobs to be run in parallel.
#     snp_type, vcffile, cur_ref, globs = select_item;
#     if snp_type == "SNP":
#         snp_vcffile, exit_flag = selectSNPs(snp_type, vcffile, cur_ref, globs);
#     elif globs['iteration'] == globs['num-iters'] and snp_type == "INDEL":
#         snp_vcffile, exit_flag = selectSNPs(snp_type, vcffile, cur_ref, globs);
#     else:
#         snp_vcffile, exit_flag = False, False;

#     return snp_type, snp_vcffile, exit_flag;

# #############################################################################

# def selectSNPs(snp_type, vcffile, cur_ref, globs):
# # Run the command to select only SNPs from a VCF file.
#     if snp_type == "SNP":
#         cur_logfile = os.path.join(globs['logdir'], "gatk-selectvariants-snps-iter-" + globs['iter-str'] + ".log");
#         snp_vcffile = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-snps.vcf.gz");
#     elif snp_type == "INDEL":
#         cur_logfile = os.path.join(globs['logdir'], "gatk-selectvariants-indels-iter-" + globs['iter-str'] + ".log");
#         snp_vcffile = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-indels.vcf.gz");

#     run_flag = PC.runCheck([snp_vcffile], cur_logfile, globs);

#     if run_flag:
#         gatk_cmd = globs['gatk-path'] + " SelectVariants -R " + cur_ref + " -V " + vcffile + " -O " + snp_vcffile + " -select-type " + snp_type;
#         exit_flag = PC.runCMD(gatk_cmd, "GATK SelectVariants", cur_logfile, True, globs);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# VCF file already exists", globs['pad'], sep=".") + snp_vcffile + "\n");
#         exit_flag = False;

#     return snp_vcffile, exit_flag;


# def varFilter(filter_item):
# # Run the command to filter variants from a VCF file based on input filters. Default: "MQ < 30.0 || DP < 5 || DP > 60"
#     snp_type, vcffile, cur_ref, globs = filter_item
    
#     if snp_type == "SNP":
#         cur_logfile = os.path.join(globs['logdir'], "gatk-varfilter-snps-iter-" + globs['iter-str'] + ".log");
#         filter_vcffile = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-snps-filtered.vcf.gz");
#     elif snp_type == "INDEL":
#         cur_logfile = os.path.join(globs['logdir'], "gatk-varfilter-indels-iter-" + globs['iter-str'] + ".log");
#         filter_vcffile = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-indels-filtered.vcf.gz");        

#     run_flag = PC.runCheck([filter_vcffile], cur_logfile, globs);

#     if run_flag:
#         gatk_cmd = globs['gatk-path'] + " VariantFiltration -R " + cur_ref + " -V " + vcffile + " -filter " + globs['filter'] + " -filter-name 'pseudoit-filter' -O " + filter_vcffile;
#         exit_flag = PC.runCMD(gatk_cmd, "GATK VariantFiltration", cur_logfile, True, globs);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# VCF file already exists", globs['pad'], sep=".") + filter_vcffile + "\n");
#         exit_flag = False;

#     return snp_type, filter_vcffile, exit_flag;

