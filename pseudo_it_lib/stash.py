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


# #############################################################################

# def genConsensus(vcffile, cur_ref, globs):
# # Run the command to generate a consensus FASTA file from the reference and the variants.
#     cur_logfile = os.path.join(globs['logdir'], "gatk-genconsensus-iter-" + globs['iter-str'] + ".log");
#     consensus_file = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + ".fa");

#     run_flag = PC.runCheck([consensus_file], cur_logfile, globs);

#     if run_flag:
#         gatk_cmd = globs['gatk-path'] + " FastaAlternateReferenceMaker -R " + cur_ref + " -V " + vcffile + " -O " + consensus_file;
#         exit_flag = PC.runCMD(gatk_cmd, "GATK FastaAlternateReferenceMaker", cur_logfile, True, globs);
#         #bcftools_cmd = globs['bcftools-path'] + " consensus -f " + cur_ref + " -o " + consensus_file + " -e \"FILTER='pseudoit-filter'\" " + vcffile
#         #exit_flag = PC.runCMD(bcftools_cmd, "bcftools consensus", cur_logfile, True, globs);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Consensus FASTA already exists", globs['pad'], sep=".") + consensus_file + "\n");
#         exit_flag = False;

#     return consensus_file, exit_flag;

# #############################################################################

# def fixHeaders(fafile, globs):
# # Fix the headers from the consensus FASTA file.
#     fa_fixed = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-final.fa");
#     sed_file = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + "-fixheaders.sed");
#     exit_flag = False;

#     if not os.path.isfile(fa_fixed) or globs['overwrite']:
#         # cur_scaffs, null = piref.getScaffs(fafile, globs, None, False);
#         # scaff_rep = {};
#         # for scaff in cur_scaffs:
#         #     scaff_orig = scaff[scaff.index(" ")+1:scaff.rindex(":")];
#         #     scaff_rep[scaff] = scaff_orig;


#         # with open(sed_file, "w") as sf:
#         #     for scaff in scaff_rep:
#         #         sf.write("s/" + scaff + "/" + scaff_rep[scaff] + "/g;\n");

#         # cmd = "sed -f " + sed_file + " < " + fafile + " > " + fa_fixed;

#         cmd = "sed 's/>[0-9]* />/g; s/:1-[0-9]*//g' " + fafile + " > " + fa_fixed;

#         if globs['dryrun']:
#             PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : " + cmd);
#         else:
#             PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : " + cmd);
#             os.system(cmd);

#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Fixed FASTA already exists", globs['pad'], sep=".") + fa_fixed + "\n");

#     return fa_fixed, exit_flag;



# #############################################################################

# def maskFa(maskbed, cur_ref, globs):
# # Fix the headers from the consensus FASTA file.
#     cur_logfile = os.path.join(globs['logdir'], "bedtools-maskfasta-" + globs['iter-str'] + ".log");
#     mask_fa = os.path.join(globs['iterdir'], "iter-" + str(int(globs['iter-str']) - 1) + "-masked.fa");

#     run_flag = PC.runCheck([mask_fa], cur_logfile, globs);

#     if run_flag:
#         mask_cmd = globs['bedtools-path'] + " maskfasta -fi " + cur_ref + " -bed " + maskbed + " -soft -fo " + mask_fa;
#         exit_flag = PC.runCMD(mask_cmd, "tabix", cur_logfile, True, globs);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Masked FASTA file already exists", globs['pad'], sep=".") + mask_fa + "\n");
#         exit_flag = False;

#     return mask_fa, exit_flag;

# #############################################################################

# def rmNonRef(vcffile, globs):
# # Deals with the weird <NON_REF> alleles in gatk's new GVCF format.
#     vcffile_decom = vcffile.replace(".gz", "");

#     cmd = "gunzip " + vcffile;
#     if globs['dryrun']:
#         PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : " + cmd);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : " + cmd);
#         os.system(cmd);

#     cmd = "sed -i 's/,<NON_REF>//g;s/<NON_REF>/./g;' " + vcffile_decom;
#     if globs['dryrun']:
#         PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : " + cmd);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : " + cmd);
#         os.system(cmd);

#     cmd = "bgzip " + vcffile_decom;
#     if globs['dryrun']:
#         PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : " + cmd);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : " + cmd);
#         os.system(cmd);

# #############################################################################

# def finalConsensus(vcffile, cur_ref, globs):
# # Run the command to generate a consensus FASTA file from the reference and the variants.
#     if globs['indels']:
#         suffix = "-indels";
#     else:
#         suffix = "-noindels";

#     cur_logfile = os.path.join(globs['logdir'], "bcftools-consensus-iter-" + globs['iter-str'] + suffix + ".log");
#     chain_file = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + suffix + "-final.chain");
#     consensus_file = os.path.join(globs['iterdir'], "iter-" + globs['iter-str'] + suffix + "-final.fa");

#     if not globs['dryrun']:
#         first_lower, linestr_orig, linestr_repl = getConsCase(cur_ref);
#         # This first_lower stuff is a hack to deal with bcftools consensus using the case of the first base in the reference fasta to inject variants.
#         # Possibly resolved: https://github.com/samtools/bcftools/issues/1150#issuecomment-582407490
#         # Need to test and make sure it is in official release before I remove this hack.
#     else:
#         first_lower, linestr_orig, linestr_repl = True, "a", "A";

#     if first_lower:
#         cmd = "sed -i '2 s/" + linestr_orig + "/" + linestr_repl + "/g' " + cur_ref;
#         if globs['dryrun']:
#             PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : sed to convert case of first ref char." );
#         else:
#             PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : sed to convert case of first ref char." );
#             os.system(cmd);
#     # Part of first_lower hack.

#     run_flag = PC.runCheck([consensus_file], cur_logfile, globs);

#     if run_flag:
#         bcftools_cmd = globs['bcftools-path'] + " consensus -f " + cur_ref + " -o " + consensus_file
#         if globs['indels']:
#             bcftools_cmd += " -c " + chain_file
#         bcftools_cmd += " -e \"FILTER='pseudoit'\" " + vcffile;
#         exit_flag = PC.runCMD(bcftools_cmd, "bcftools consensus", cur_logfile, True, globs);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Consensus FASTA already exists", globs['pad'], sep=".") + consensus_file + "\n");
#         exit_flag = False;
#     # Consensus command

#     if first_lower:
#         if not globs['dryrun']:
#             first_lower, linestr_orig, linestr_repl = getConsCase(cur_ref);

#         cmd = "sed -i '2 s/" + linestr_orig + "/" + linestr_repl + "/g' " + cur_ref;
#         if globs['dryrun']:
#             PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : sed to revert case of first consensus char.");
#         else:
#             PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : sed to revert case of first consensus char.");
#             os.system(cmd);        

#         if not globs['dryrun']:
#             first_lower, linestr_orig, linestr_repl = getConsCase(consensus_file);

#         cmd = "sed -i '2 s/" + linestr_orig + "/" + linestr_repl + "/g' " + consensus_file;
#         if globs['dryrun']:
#             PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Would execute   : sed to revert case of first ref char.");
#         else:
#             PC.printWrite(globs['logfilename'], globs['log-v'], "# " + PC.getDateTime() + " --> Executing   : sed to revert case of first ref char.");
#             os.system(cmd);
#     # Part of first_lower hack.

#     return consensus_file, exit_flag;

# #############################################################################

# def getConsCase(cur_ref):
# # bcftools consensus uses the case of the first base in the reference file for injecting
# # variants. This is a hack to make sure it is always upper case.

#     first_lower = False;
#     lines = 0;
#     for line in open(cur_ref):
#         line = line.strip();
#         lines += 1;
#         if lines == 2:
#             linestr_orig = line;
#             linelist = list(linestr_orig);
#             if linelist[0].islower():
#                 linelist[0] = linelist[0].upper();
#                 first_lower = True;
#             elif linelist[0].isupper():
#                 linelist[0] = linelist[0].lower();
#             linestr_repl = "".join(linelist);
#             break;

#     return first_lower, linestr_orig, linestr_repl;

# #############################################################################


# def report_stats(globs, msg="", step_start=0, stat_start=False, stat_end=False, sep=" "):
# Uses psutil to gather memory and time info between steps and print them to the screen.
	# import timeit
	# if globs['psutil']:
	# 	import psutil;
	# 	dashes = 161;
	# else:
	# 	dashes = 125;
	# cur_time = timeit.default_timer();
	# if stat_start:
	# # The first time through just print the headers.
	# 	globs['progstarttime'] = cur_time;
	# 	printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * dashes);
	# 	if globs['psutil']:
	# 		printWrite(globs['logfilename'], globs['log-v'], "# Date/time" + " " * 13 + "Current step" + " " * 25 + "Time since prev (sec)" + " " * 6 + "Elapsed time (sec)" + " " * 4 + "Current mem usage (MB)" + " " * 4 + "Virtual mem usage (MB)");
	# 	else:
	# 		printWrite(globs['logfilename'], globs['log-v'], "# Date/time" + " " * 13 + "Current step" + " " * 25 + "Time since prev (sec)" + " " * 6 + "Elapsed time (sec)");
	# 	printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * dashes);
	# else:
	# 	prog_elapsed = round(cur_time - globs['progstarttime'], 5);
	# 	step_elapsed = round(cur_time - step_start, 5);
	# 	if globs['psutil']:
	# 		mem = round(sum([p.memory_info()[0] for p in globs['pids']]) / float(2 ** 20), 5);
	# 		vmem = round(sum([p.memory_info()[1] for p in globs['pids']]) / float(2 ** 20), 5);
	# 		printWrite(globs['logfilename'], globs['log-v'], "# " + getDateTime() + " " + msg + sep * (37-len(msg)) + str(step_elapsed) + sep * (27-len(str(step_elapsed))) + str(prog_elapsed) + sep * (22-len(str(prog_elapsed))) + str(mem) + sep * (26-len(str(mem))) + str(vmem));
	# 	else:
	# 		printWrite(globs['logfilename'], globs['log-v'], "# " + getDateTime() + " " + msg + sep * (37-len(msg)) + str(step_elapsed) + sep * (27-len(str(step_elapsed))) + str(prog_elapsed) + sep * (22-len(str(prog_elapsed))));
	# 	if stat_end:
	# 		printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * dashes);
	# return cur_time;



# def indexVCF(vcffile, globs, suffix=""):
# # Index the combined VCF from gatherVcfs.

#     cur_logfile = os.path.join(globs['iterlogdir'], "vcf-index-iter-" + globs['iter-str'] + ".log");

#     if suffix != "":
#         suffix = "-" + suffix;
#     cur_logfile = os.path.join(globs['logdir'], "vcf-index-iter-" + globs['iter-str'] + suffix + ".log");
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

# #############################################################################

# def getSubPID(n):
# # Gets the process ids for the --stats option.
# 	import psutil
# 	return psutil.Process(os.getpid());



# def indexDistributor(index_item):
# # This function distributes the FASTA indexing jobs to be run in parallel.
#     cur_opt, cur_ref, globs, step_start_time = index_item;
#     if cur_opt == 'dict':
#         step_start_time, exit_flag = picardDict(cur_ref, globs, step_start_time);
#     elif cur_opt == 'faidx':
#         step_start_time, exit_flag = faidxFasta(cur_ref, globs, step_start_time);
#     elif cur_opt == 'index':
#         step_start_time, exit_flag = indexFasta(cur_ref, globs, step_start_time);
    
#     return step_start_time, exit_flag;

# ############################################################################

# def picardDict(cur_ref, globs, step_start_time):
# # Create the sequence dictionary for the FASTA.    
#     ref_ext = PC.detectRefExt(cur_ref, globs);
#     cur_logfile = os.path.join(globs['logdir'], "picard-dict-iter-" + globs['iter-str'] + ".log");
#     dictfile = cur_ref.replace(ref_ext, ".dict");

#     run_flag = PC.runCheck([dictfile], cur_logfile, globs);

#     if run_flag:
#         #step_start_time = PC.report_stats(globs, "--> ITERATION " + globs['iter-str'] + ": Dictionary", step_start=step_start_time);
#         dict_cmd = globs['picard-path'] + " CreateSequenceDictionary R=" + cur_ref + " O=" + dictfile;
#         exit_flag = PC.runCMD(dict_cmd, "Picard CreateSequenceDictionary", cur_logfile, True, globs);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# Dictionary file already exists", globs['pad'], sep=".") + dictfile);
#         exit_flag = False;
#     return step_start_time, exit_flag;

# #############################################################################

# def faidxFasta(cur_ref, globs, step_start_time):
# # Index the FASTA with samtools faidx.
#     cur_logfile = os.path.join(globs['logdir'], "samtools-faidx-iter-" + globs['iter-str'] + ".log");
#     faidxfile = cur_ref + ".fai";

#     run_flag = PC.runCheck([faidxfile], cur_logfile, globs);

#     if run_flag:        
#         #step_start_time = PC.report_stats(globs, "--> ITERATION " + globs['iter-str'] + ": faidx", step_start=step_start_time);
#         faidx_cmd = globs['samtools-path'] + " faidx " + cur_ref;
#         exit_flag = PC.runCMD(faidx_cmd, "samtools faidx", cur_logfile, True, globs);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# faidx file already exists", globs['pad'], sep=".") + faidxfile);
#         exit_flag = False;
#     # samtools faidx
#     return step_start_time, exit_flag;

# #############################################################################

# def indexFasta(cur_ref, globs, step_start_time):
# # Index the FASTA file with BWA.
#     cur_logfile = os.path.join(globs['logdir'], "bwa-index-iter-" + globs['iter-str'] + ".log");
#     indexfiles = [cur_ref + ".amb", cur_ref + ".ann", cur_ref + ".bwt", cur_ref + ".pac", cur_ref + ".sa"];

#     run_flag = PC.runCheck(indexfiles, cur_logfile, globs);

#     if run_flag:
#         #step_start_time = PC.report_stats(globs, "--> ITERATION " + globs['iter-str'] + ": BWA index", step_start=step_start_time);
#         index_cmd = globs['bwa-path'] + " index " + cur_ref;
#         exit_flag = PC.runCMD(index_cmd, "BWA index", cur_logfile, True, globs);
#     else:
#         PC.printWrite(globs['logfilename'], globs['log-v'], PC.spacedOut("# BWA index files already exist", globs['pad'], sep=".") + ",".join(indexfiles));
#         exit_flag = False;
#     # bwa index
#     return step_start_time, exit_flag;

#############################################################################

	# if args.java_heap and not isdigit(args.java_heap) or int(args.java_heap) < 1:
	# 	PC.errorOut("OP8", "-heap must be a positive integer.", globs);
	# elif args.java_heap:
	# 	globs['heap'] = int(args.java_heap);
	# Getting the java heap option.


    	# step_start_time = "";
	# if globs['psutil']:
	# 	globs['pids'] = [psutil.Process(os.getpid())];	
	# globs['stats'] = True;
	# if not globs['norun']:
	# 	step_start_time = PC.report_stats(globs, stat_start=True);
	# Initializing the stats options if --quiet is not set.

#############################################################################

# def selectSNPsScaff(globs, cmds):
# # Run the command to select only SNPs from a VCF file.

#     gatk_cmds = {};
#     for scaff in globs['scaffolds']:
#         cur_logfile = os.path.join(globs['itervcflogdir'], "gatk-selectvariants-iter-" + scaff + "-" + globs['iter-str'] + ".log")
#         vcf_file = os.path.join(globs['itervcfscaffdir'], scaff + "-iter-" + globs['iter-str'] + ".vcf.gz");
#         snp_file = os.path.join(globs['itervcfscaffdir'], scaff + "-iter-" + globs['iter-str'] + "-snps.vcf.gz");

#         gatk_cmd = globs['gatk-path'] + " SelectVariants -V " + vcf_file + " -O " + snp_file + " -select-type SNP -xl-select-type INDEL -xl-select-type MIXED -xl-select-type SYMBOLIC";

#         cmd_num = PC.getCMDNum(globs, len(cmds));
#         cmds[gatk_cmd] = { 'cmd-num' : cmd_num, 'desc' : "Select SNPs " + scaff, 'outfile' : vcf_file,  'logfile' : cur_logfile, 'start' : False };
#         gatk_cmds[gatk_cmd] = { 'cmd-num' : cmd_num, 'desc' : "Select SNPs " + scaff, 'outfile' : vcf_file,  'logfile' : cur_logfile, 'start' : False };

#     if globs['dryrun']:
#         cmd_num = PC.getCMDNum(globs, len(cmds));
#         gatk_skeleton_cmd = globs['gatk-path'] + " SelectVariants -R <reference FASTA> -V <input vcf> -O <snp only vcf> -select-type SNP -xl-select-type INDEL -xl-select-type MIXED -xl-select-type SYMBOLIC";
#         cmds[gatk_skeleton_cmd] = { 'cmd-num' : cmd_num, 'desc' : str(globs['num-procs']) + " SelectVariants procs in parallel", 'outfile' : "",  'logfile' : "", 'start' : False };
#         PC.report_step(globs, cmds, gatk_skeleton_cmd, "DRYRUN", gatk_skeleton_cmd);

#     else:
#         pool = mp.Pool(processes=globs['num-procs']);
#         for result in pool.starmap(PC.runCMD, ((gatk_cmd, globs, cmds, True) for gatk_cmd in gatk_cmds )):
#             if result:
#                 pool.terminate();
#                 globs['exit-code'] = 1;
#                 PC.endProg(globs);
#         pool.terminate();    

#     return cmds;

#############################################################################

	# try:
	# 	import argparse;
	# except:
	# 	PC.errorOut("\n*** ERROR: Your installation of Python is missing the argparse module. Please try a different version of Python (3+) or install the module.\n")
	# # First check if the argparse module is installed. If not, the input options cannot be parsed.