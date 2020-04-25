# Iterative mapping steps for Pseudo-it
#############################################################################
import sys, os, subprocess, lib.picore as PC, lib.piref as piref, \
    lib.pimap as pimap, lib.varprep as varprep, multiprocessing as mp, \
    lib.varcall as varcall, lib.varpost as varpost
#############################################################################

def mapping(globs, step_start_time):

    globs = PC.getIterStr(globs);
    pool = mp.Pool(processes=globs['num-procs']);
    cur_ref = PC.getRef(globs);
    # Iteration prep

    globs['iterdir'] = os.path.join(globs['outdir'], "iter-" + globs['iter-str']);
    if not os.path.isdir(globs['iterdir']):
        os.makedirs(globs['iterdir']);
    # Make directory for current iteration.

    if globs['iteration'] == 1:
        globs['scaffolds'], step_start_time = piref.getScaffs(globs['ref'], globs, step_start_time);

        piref.indexCheck(globs['ref'], globs);
    # Get the scaffold IDs from the reference FASTA and check that all index files have been created.

    if globs['iteration'] != 1:
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Prep FASTA", step_start=step_start_time);
        for result in pool.imap(piref.indexDistributor, ((index_opt, cur_ref, globs, step_start_time) for index_opt in ['dict', 'faidx', 'index'])):
            step_start_time, exit_flag = result;
            if exit_flag:
                pool.terminate();
                PC.endProg(globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    # Prepare the current reference FASTA -- indexes with samtools, bwa, and creates picard dictionary.

    step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Mapping", step_start=step_start_time);
    cur_bamfiles = [];
    for result in pool.imap(pimap.BWA, ((cur_lib, globs['libs'][cur_lib], cur_ref, globs, step_start_time) for cur_lib in globs['libs'])):
        exit_flag = result[1];
        if exit_flag:
            pool.terminate();
            PC.endProg(globs);
        cur_bamfiles.append(result[0]);
    PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    # Map reads in parallel for all given libraries.

    if len(cur_bamfiles) > 1:
        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": BAM Merging", step_start=step_start_time);    
        merged_bamfile, exit_flag = pimap.mergeBam(cur_bamfiles, globs);
        PC.exitCheck(exit_flag, globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    else:
        merged_bamfile = cur_bamfiles[0];
    # Merge BAM files if more than one library.
   

    step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Add Read groups", step_start=step_start_time);
    rg_bamfile, exit_flag = varprep.addRG(merged_bamfile, globs);
    PC.exitCheck(exit_flag, globs);
    PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    # Add read groups.

    step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Mark duplicates", step_start=step_start_time);
    dup_bamfile, exit_flag = varprep.markDups(rg_bamfile, globs);
    PC.exitCheck(exit_flag, globs);
    PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    # Mark duplicates.

    step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Index BAM", step_start=step_start_time);
    exit_flag = varprep.indexBAM(dup_bamfile, globs);
    PC.exitCheck(exit_flag, globs);
    PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    # Index BAM file.

    step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": HaplotypeCaller", step_start=step_start_time);
    if globs['num-procs'] == 1 or len(globs['scaffolds']) == 1:
        vcf_file, exit_flag = varcall.haplotypeCaller(dup_bamfile, cur_ref, globs);
        PC.exitCheck(exit_flag, globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # HaplotypeCaller.

        if globs['iteration'] == globs['num-iters']:
            lastIter(vcf_file, "", cur_ref, globs, step_start_time);
        # If this is the last iteration, go to those steps.
    # Serial version.
    ###
    else:
        PC.printWrite(globs['logfilename'], globs['log-v'], "# --> Running " + str(globs['num-procs']) + " scaffolds in parallel.");
        vcf_dir = os.path.join(globs['iterdir'], "vcf-scaffolds");
        if not os.path.isdir(vcf_dir):
            os.makedirs(vcf_dir);
        vcf_log_dir = os.path.join(globs['iterdir'], "vcf-scaffolds-logs");
        if not os.path.isdir(vcf_log_dir):
            os.makedirs(vcf_log_dir);
        # Make directories for the VCFs.

        for result in pool.imap(varcall.haplotypeCallerMulti, ((scaff, dup_bamfile, cur_ref, vcf_dir, vcf_log_dir, globs, step_start_time) for scaff in globs['scaffolds'])):
            step_start_time, exit_flag = result;
            if exit_flag:
                pool.terminate();
                PC.endProg(globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # HaplotypeCaller.

        if globs['iteration'] == globs['num-iters']:
            pool.terminate();
            lastIter(vcf_dir, vcf_log_dir, cur_ref, globs, step_start_time);
        # If this is the last iteration, go to those steps.

        else:
            step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": CatVariants", step_start=step_start_time);
            vcf_file, exit_flag = varcall.gatherVcfs(vcf_dir, cur_ref, globs);
            PC.exitCheck(exit_flag, globs);
            PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
            # Combine the scaffold VCF files.

            step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Index VCF", step_start=step_start_time);
            exit_flag = varcall.indexVCF(vcf_file, globs);
            PC.exitCheck(exit_flag, globs);
            PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
            # Index the combined VCF file.
    # Parallel version by scaffold.
    ###
    ## Call variants.

    if globs['iteration'] != globs['num-iters']:
    # Only do these steps if its not the last iteration.
        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Select SNPs", step_start=step_start_time);
        vcf_snp_file, exit_flag = varpost.selectSNPs(vcf_file, cur_ref, globs);
        PC.exitCheck(exit_flag, globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # Select SNPs only.

        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Filter variants", step_start=step_start_time);
        vcf_filter_file, exit_flag = varpost.varFilter((vcf_snp_file, "", "", cur_ref, globs));
        PC.exitCheck(exit_flag, globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");    
        # Filter VCF

        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Index filtered VCF", step_start=step_start_time);
        exit_flag = varcall.indexVCF(vcf_filter_file, globs, suffix="filtered");
        PC.exitCheck(exit_flag, globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # Index the filtered VCF file.

        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Generate consensus", step_start=step_start_time);
        cur_fa, exit_flag = varpost.genConsensus(vcf_filter_file, cur_ref, globs);
        PC.exitCheck(exit_flag, globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # Generate conesnsus sequence between refrence and variant calls

        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Fix headers", step_start=step_start_time);
        cur_fa = varpost.fixHeaders(cur_fa, globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # Put back the original FASTA headers in the conensus

        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + " Complete!", step_start=step_start_time);
        PC.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 100);
        pool.terminate();
        # Finish up the iteration.

    globs['iteration'] +=1;
    return globs, step_start_time;

#############################################################################

def lastIter(gvcf, gvcf_log, cur_ref, globs, step_start_time):

    pool = mp.Pool(processes=globs['num-procs']);
    step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": GenotypeGVCFs", step_start=step_start_time);
    if globs['num-procs'] == 1 or len(globs['scaffolds']) == 1:
        vcf_file, exit_flag = varcall.genotypeGVCFs(gvcf, cur_ref, globs);
        PC.exitCheck(exit_flag, globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # Genotype GVCFs.

        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Filter variants", step_start=step_start_time);
        vcf_file, exit_flag = varpost.varFilter((vcf_file, "", "", cur_ref, globs));
        PC.exitCheck(exit_flag, globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # Filter variants.    
    # Serial version.
    ###
    else:
        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": GenotypeGVCFs", step_start=step_start_time);
        PC.printWrite(globs['logfilename'], globs['log-v'], "# --> Running " + str(globs['num-procs']) + " scaffolds in parallel.");
        for result in pool.imap(varcall.genotypeGVCFsMulti, ((scaff, cur_ref, gvcf, gvcf_log, globs, step_start_time) for scaff in globs['scaffolds'])):
            step_start_time, exit_flag = result;
            if exit_flag:
                pool.terminate();
                PC.endProg(globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # Genotype GVCFs.

        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Filter variants", step_start=step_start_time);
        PC.printWrite(globs['logfilename'], globs['log-v'], "# --> Running " + str(globs['num-procs']) + " scaffolds in parallel.");
        vcf = [ os.path.join(gvcf, f) for f in os.listdir(gvcf) if f.endswith(".vcf.gz")];
        for result in pool.imap(varpost.varFilter, ((scaff_vcf, gvcf, gvcf_log, cur_ref, globs) for scaff_vcf in vcf)):
            vcf_filter_file, exit_flag = result;
            if exit_flag:
                pool.terminate();
                PC.endProg(globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # Filter variants.       

        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": GatherVcfs", step_start=step_start_time);
        vcf_file, exit_flag = varcall.gatherVcfs(gvcf, cur_ref, globs);
        PC.exitCheck(exit_flag, globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # Combine filtered VCFs.
    # Parallel version by scaffold.
    ###
    ## Genotype variants.

    # step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Remove <NON_REF>", step_start=step_start_time);
    # varpost.rmNonRef(vcf_file, globs);
    # PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    # # Because of a bug in GenotypeGVCFs, some of the weird <NON_REF> alleles remain in the final VCF. This messes up bcftools consensus.
    # # This hack removes the <NON_REF>s.
    # # Supposedly this will be fixed in the next GATK release:
    # # https://gatk.broadinstitute.org/hc/en-us/community/posts/360056352871-Some-NON-REF-alleles-remain-after-GenotypeGVCFs-when-using-include-non-variant-sites
    ### They seem to have fixed this in version 4.1.5.0

    step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Index filtered VCF", step_start=step_start_time);
    exit_flag = varcall.indexVCF(vcf_file, globs, suffix="filtered");
    PC.exitCheck(exit_flag, globs);
    PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    # Index filtered VCF.

    step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Get mask sites", step_start=step_start_time);
    mask_bed, exit_flag = varpost.getMask(vcf_file, globs);
    PC.exitCheck(exit_flag, globs);
    PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    # Get mask sites from final VCF.

    step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Mask previous iteration FASTA", step_start=step_start_time);
    mask_fa, exit_flag = varpost.maskFa(mask_bed, cur_ref, globs);
    PC.exitCheck(exit_flag, globs);
    PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    # Mask the current (last iteration) reference FASTA.

    if not globs['indels']:
        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Select SNPs", step_start=step_start_time);
        vcf_file, exit_flag = varpost.selectSNPs(vcf_file, cur_ref, globs);
        PC.exitCheck(exit_flag, globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # Select SNPs only.

        step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Index filtered SNP VCF", step_start=step_start_time);
        exit_flag = varcall.indexVCF(vcf_file, globs, suffix="filtered-snps");
        PC.exitCheck(exit_flag, globs);
        PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
        # Index the filtered VCF SNP file.

    step_start_time = PC.report_stats(globs, "ITERATION " + globs['iter-str'] + ": Generate final FASTA", step_start=step_start_time);
    mask_fa, exit_flag = varpost.finalConsensus(vcf_file, mask_fa, globs);
    PC.exitCheck(exit_flag, globs);
    PC.printWrite(globs['logfilename'], globs['log-v'], "\n");
    # Make the final consensus FASTA from the masked FASTA and the final filtered VCF file.

    pool.terminate();