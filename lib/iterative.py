# Iterative mapping steps for Pseudo-it
#############################################################################
import sys, os, shutil, subprocess, lib.picore as PC, lib.piref as piref, \
    lib.pimap as pimap, multiprocessing as mp, \
    lib.varcall as varcall, lib.varpost as varpost, lib.consensus as con
#############################################################################

def mapping(globs):

    globs = PC.getIterStr(globs);
    PC.printWrite(globs['logfilename'], globs['log-v'], "#\n# " + "=" * 51 + " ITERATION " + globs['iter-str'] + " STARTING! " + "=" * 50);
    # Update the current iteration number in the global dict and print a status report.

    globs['iterstarttime'] = PC.report_step(globs, "", start=True);
    if globs['iteration'] == 1:
        globs['progstarttime'] = globs['iterstarttime'];
    # If it is the first iteration, save the start time as the overall start time.

    cmds = {};
    # A dictionary of all commands that are run throughout the iteration.
    
    cur_ref = PC.getRef(globs);
    if globs['iteration'] == globs['num-iters']:
        globs['last-iter'] = True;
    # Iteration prep: Get the correct reference file name and check if we're on the last iteration.

    ##########

    globs['iterdir'] = os.path.join(globs['outdir'], "iter-" + globs['iter-str']);
    globs['iterbamdir'] = os.path.join(globs['iterdir'], "bam");
    globs['itervcfdir'] = os.path.join(globs['iterdir'], "vcf");
    globs['iterfadir'] = os.path.join(globs['iterdir'], "fa");
    globs['iterlogdir'] = os.path.join(globs['iterdir'], "logs");
    for d in [ globs['iterdir'], globs['iterbamdir'], globs['itervcfdir'], globs['iterfadir'], globs['iterlogdir'] ]:
        if not os.path.isdir(d):
            if globs['map-only'] and d not in [ globs['iterbamdir'], globs['iterlogdir'] ]:
                continue;
            os.makedirs(d);
    # Each iteration consists of an overall directory (iter-N), and sub-directories for bam, vcf, fa, and logs.

    ##########

    if globs['in-bed'] and globs['bed-mode'] == "file":
        globs['itervcfscaffdir'] = globs['itervcfdir'];
        globs['itervcflogdir'] = globs['itervcfdir'];  
    # If a bed file has been provided in with -bedmode file, no further subdirectories are needed.    
    elif globs['last-iter'] and globs['mask'] != "none":
        globs['itervcfscaffdir'] = os.path.join(globs['itervcfdir'], "gvcf-scaff");
        globs['itervcflogdir'] = os.path.join(globs['itervcfdir'], "gvcf-logs");
    else:
        globs['itervcfscaffdir'] = os.path.join(globs['itervcfdir'], "vcf-scaff");
        globs['itervcflogdir'] = os.path.join(globs['itervcfdir'], "vcf-logs");
    if not globs['map-only'] and globs['bed-mode'] in [False, "regions"]:
        for d in [ globs['itervcfscaffdir'], globs['itervcflogdir'] ]:
            if not os.path.isdir(d):
                os.makedirs(d);
    # Without a bed file in -bedmode file,  the vcf directory for the current iteration also needs a subdirectory for the separate vcf files for
    # each interval to run in parallel. For the last iteration these are gvcfs if masking.

    ##########

    if globs['bam'] and globs['iteration'] == 1:
        globs['iter-final-bam-log'] = "NA";
        globs['iter-final-bam'] = globs['bam'];  
    # If a bam file is provided with -bam, set that as the final expected bam file for the first iteration.
    else:
        if globs['mkdups']:
            globs['iter-final-bam-log'] = os.path.join(globs['iterlogdir'], "picard-mkdup-iter-" + globs['iter-str'] + ".log");
            globs['iter-final-bam'] = os.path.join(globs['iterbamdir'], "merged-rg-mkdup-iter-" + globs['iter-str'] + ".bam.gz");
        else:
            globs['iter-final-bam-log'] = os.path.join(globs['iterlogdir'], "picard-merge-bam-iter-" + globs['iter-str'] + ".log");
            globs['iter-final-bam'] = os.path.join(globs['iterbamdir'], "merged-iter-" + globs['iter-str'] + ".bam.gz");
        # If --nomkdup is set, the final BAM file for each iteration should not have the mkdup suffix.        
    # The expected final BAM file for this iteration.

    ##########

    if globs['last-iter']:
        globs['iter-gather-vcf-log'] = os.path.join(globs['iterlogdir'], "gatk-gathervcfs-iter-" + globs['iter-str'] + ".log");
        if globs['in-bed'] and globs['bed-mode'] == "file":
            globs['iter-gather-vcf'] = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + "-filter.vcf.gz");
        else:
            globs['iter-gather-vcf'] = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + "-filter.vcf.gz");

        if globs['indels']:
            globs['iter-final-vcf-log'] = globs['iter-gather-vcf-log'];
            globs['iter-final-vcf'] = globs['iter-gather-vcf'];
        else:
            globs['iter-final-vcf-log'] = os.path.join(globs['iterlogdir'], "gatk-selectsnps-iter-" + globs['iter-str'] + ".log");
            globs['iter-final-vcf'] = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + "-filter-snps.vcf.gz");
    else:
        globs['iter-gather-vcf-log'] = os.path.join(globs['iterlogdir'], "gatk-gathervcfs-iter-" + globs['iter-str'] + "-intermediate.log");
        if globs['in-bed'] and globs['bed-mode'] == "file":
            globs['iter-gather-vcf'] = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + "-filter.vcf.gz");
        else:
            globs['iter-gather-vcf'] = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + "-filter-intermediate.vcf.gz");

        globs['iter-final-vcf-log'] = os.path.join(globs['iterlogdir'], "gatk-selectsnps-iter-" + globs['iter-str'] + "-intermediate.log");
        globs['iter-final-vcf'] = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + "-filter-intermediate-snps.vcf.gz");
    globs['intermediate-vcf'] = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + ".vcf.gz");
    globs['intermediate-vcf-log'] = os.path.join(globs['itervcfdir'], "gatk-haplotypcaller-iter-" + globs['iter-str'] + ".log");
    # The expected final vcf file for each iteration. This differs depending on if it is the last iteration, whether a bed
    # file is provided with -bedmode file, and whether --noindels is specified.

    ##########

    if globs['last-iter']:
        fa_suffix = "";
        if not globs['indels']:
            fa_suffix += "-snps";
        if globs['diploid']:
            fa_suffix += "-diploid";
        if globs['mask'] in ["soft", "hard"]:
            fa_suffix += "-" + globs['mask'] + "mask";
        else:
            fa_suffix += "-nomask";
        fa_suffix += "-final";

        globs['iter-consensus-log'] = os.path.join(globs['iterlogdir'], "bcftools-consensus-iter-" + globs['iter-str'] + fa_suffix + ".log");
        globs['iter-final-chain'] = os.path.join(globs['iterfadir'], "iter-" + globs['iter-str'] + fa_suffix +  ".chain");
        globs['iter-final-fa'] = os.path.join(globs['iterfadir'], "iter-" + globs['iter-str'] + fa_suffix +  ".fa");

    else:
        globs['iter-consensus-log'] = os.path.join(globs['iterlogdir'], "bcftools-consensus-iter-" + globs['iter-str'] + "-snps-intermediate.log");
        globs['iter-final-chain'] = os.path.join(globs['iterfadir'], "iter-" + globs['iter-str'] + "-snps-intermediate.chain");
        globs['iter-final-fa'] = os.path.join(globs['iterfadir'], "iter-" + globs['iter-str'] + "-snps-intermediate.fa");
    # The expected final FASTA files for this iteration. This depends on whether --noindels, --diploid, and -mask are specified.

    # Output files for current iteration
    ##########

    if globs['iteration'] == 1:
        cmds = piref.indexCheck(globs['ref'], globs, cmds);
        globs['scaffolds'], cmds = piref.getScaffs(globs['ref'], globs, cmds);
    # Check that all index files have been created and get the scaffold IDs from the reference FASTA

    if globs['bam'] and globs['iteration'] == 1:
        do_mapping = False;
    elif globs['resume']:
        do_mapping = PC.prevCheck(globs['iter-final-bam'], globs['iter-final-bam-log'], globs);
    else:
        do_mapping = True;
    # Determine whether mapping needs to be done for this iteration. Depends on whether -resume or -bam are specified.

    if do_mapping:
        do_varcalling = True;
    elif globs['resume']:
        do_varcalling = PC.prevCheck(globs['iter-final-vcf'], globs['iter-final-vcf-log'], globs);
    else:
        do_varcalling = True;
    # Determine whether variant calling needs to be done for this iteration. Depends on whether mapping needs to be done
    # or whether -resume is specified.

    if do_varcalling:
        do_consensus = True;
    elif globs['resume']:
        do_consensus = PC.prevCheck(globs['iter-final-fa'], globs['iter-consensus-log'], globs);
    else:
        do_consensus = True;
    # Determine whether consensus generation needs to be done for this iteration. Depends on whether variant calling needs to be done
    # or whether -resume is specified.   

    # CHECK WHICH STEPS WE NEED TO PERFORM
    ##########

    statstr = "EXECUTING";
    if globs['resume']:
        statstr = "RESUME";
    if globs['dryrun']:
        statstr = "DRYRUN";
    # Status string for the main step reports depending on the runtype.

    if globs['iteration'] != 1:
        cmds = piref.indexFa(globs, cmds, cur_ref);
    # If not the first iteration, index the previous iterations consensus sequence.
    else:
        if globs['in-bed'] and globs['bed-mode'] == "regions":
                globs['regions'], cmds = PC.readBed(globs['in-bed'], globs, cmds);
        # If a bed file has been provided to limit variant calling to specific regions, read it here.

        if globs['in-vcf']:
            globs['filter-sites'], cmds = PC.readVCF(globs['in-vcf']);
        # If a VCF file has been provided to filter SNPs with -vcf, read those SNPs here. This is slow, meant for a small number of SNPs
        # in small genomes only.
    # If it is the first iteration, read the supplemental files here if provided.

    ##########

    if not globs['in-bed']:
        globs['intervals'] = globs['scaffolds'];
    else:
        if globs['bed-mode'] == "regions":
            globs['intervals'] = globs['regions'];
        elif globs['bed-mode'] == "file":
            globs['intervals'] = [globs['in-bed']];
    # Retrieve the intervals over which to call variants. Depends on whether a bed file is provided and the -bedmode.

    ##########

    if do_mapping:
        PC.report_step(globs, cmds, "NA--01   Read mapping", statstr, "Mapping reads and post-processing.");

        if globs['iteration'] == 1:
            globs, cmds = pimap.getRG(globs, cmds);
        # Get the read group information if it is the first iteration.

        if globs['mapper'] == "bwa":
            bamfiles, cmds = pimap.BWA(globs, cmds, cur_ref);
        # If --mapper is bwa
        if globs['mapper'] == "hisat2":
            bamfiles, cmds = pimap.hisat2(globs, cmds, cur_ref);
        # If --mapper is hisat2
        # READ MAPPING

        #rg_bamfile, cmds = varprep.addRG(globs, cmds, bamfiles);
        # ADD READ GROUPS.
        # This is now done with the specified mapper.

        merged_bamfile, cmds = pimap.mergeBam(globs, cmds, bamfiles);
        # MERGE BAM FILES, also sorts
        
        if globs['mkdups']:
            cmds = pimap.markDups(globs, cmds, merged_bamfile);
        # MARK DUPLICATES, also generates index for merged bam file.
        
        #cmds = varprep.indexBAM(globs, cmds);
        # INDEX BAM FILE
        # Now done during MarkDuplicates or Merging
        
    elif globs['bam'] and globs['iteration'] == 1:
        PC.report_step(globs, cmds, "NA--01   Read mapping", "BAM", "initial BAM file provided, skipping all mapping steps: " + globs['iter-final-bam']);
    else:
        PC.report_step(globs, cmds, "NA--01   Read mapping", "RESUME", "previous processed BAM file found, skipping all mapping steps: " + globs['iter-final-bam']);
    
    ## READ MAPPING STEPS
    ##########

    if globs['map-only']:
        PC.report_step(globs, cmds, "NA--04   Iteration cleanup", statstr, "Removing intermediate files based on --keep* options.");
        cmds = cleanUp(globs, cmds);

        PC.printWrite(globs['logfilename'], globs['log-v'], "#\n# " + "=" * 51 + " ITERATION " + globs['iter-str'] + " COMPLETE! " + "=" * 50);
        PC.report_step(globs, "", end=True);

        globs['iteration'] +=1;
        return globs;
    # This stops the program after the first iteration of mapping if --maponly is set.

    ##########

    if do_varcalling:
        PC.report_step(globs, cmds, "NA--02   Variant calling", statstr, "Calling and post-processing variants.");

        cmds = varcall.haplotypeCaller(globs, cmds, cur_ref, globs['iter-final-bam']);
        # HAPLOTYPECALLER

        if globs['last-iter'] and globs['mask'] != "none":
            cmds = varcall.genotypeGVCFs(globs, cmds, cur_ref);
        # GENOTYPE GVCFS FOR LAST ITER

        cmds = varpost.varFilter(globs, cmds, cur_ref);
        # FILTER VCFs

        if globs['in-vcf']:
            varpost.varFilterManual(globs, cmds);
        # If a vcf file has been provided for filtering SNPs with -vcf, do that filtering here. Only good for very small genomes.

        if globs['bed-mode'] in [False, "regions"]:
            cmds = varpost.gatherVCFs(globs, cmds);
        # COMBINE VCF if a bed file hasn't been provided with -bedmode file

        cmds = varpost.indexVCF(globs, cmds, globs['iter-gather-vcf']);
        # INDEX VCF        
    else:
        PC.report_step(globs, cmds, "NA--02   Variant calling", "RESUME", "previous processed VCF file found, skipping all variant calling steps: " + globs['iter-final-vcf']);
    
    ## VARIANT CALLING STEPS
    ##########
   
    if do_consensus:
        PC.report_step(globs, cmds, "NA--03   Consensus generation", statstr, "Generating consensus FASTA.");
        if globs['last-iter'] and globs['mask'] != "none":
            mask_bedfile, cmds = con.getMask(globs, cmds, globs['iter-gather-vcf']);
            # GET MASK SITES

            cur_ref, cmds = con.maskFa(globs, cmds, mask_bedfile, cur_ref);
            # MASK PREVIOUS REFERENCE
        # Masking sites if -mask none isn't specified.

        if not globs['last-iter'] or (globs['last-iter'] and not globs['indels']):
            cmds = varpost.selectSNPs(globs, cmds, globs['iter-gather-vcf']);
            # SELECT SNPs FROM VCF IF IT IS NOT THE LAST ITERATION OR IF --noindels IS SET

            cmds = varpost.indexVCF(globs, cmds, globs['iter-final-vcf']);
            # INDEX FINAL VCF
        # For wach iteration and the last iteration if --noindels is specified, select SNPs only from the vcf here.

        cmds, globs = con.genConsensus(globs, cmds, globs['iter-final-vcf'], cur_ref);
        # GENERATE CONSENSUS
    else:
        PC.report_step(globs, cmds, "NA--03   Consensus generation", "RESUME", "previous processed consensus FASTA file found, skipping all consensus generation steps: " + globs['iter-final-fa']);    
        globs['consensus-file'] = globs['iter-final-fa'];
    
    ## CONSENSUS STEPS
    ##########

    PC.report_step(globs, cmds, "NA--04   Iteration cleanup", statstr, "Removing intermediate files based on --keep* options.");
    cmds = cleanUp(globs, cmds);

    PC.printWrite(globs['logfilename'], globs['log-v'], "#\n# " + "=" * 51 + " ITERATION " + globs['iter-str'] + " COMPLETE! " + "=" * 50);
    PC.report_step(globs, "", end=True);
    # Clean-up and status reporting

    ##########

    globs['iteration'] +=1;
    return globs;

#############################################################################

def cleanUp(globs, cmds):
    i = globs['iter-str'];
    prev_i = str(int(i) - 1);
    if len(prev_i) == 1:
        prev_iter = "0" + prev_i;

    possible_map_files = { 'iter-' + i + "-dupmets.txt" : 2,  "merged-iter-" + i + ".bam.gz" : 2, "merged-rg-iter-" + i + ".bam.gz" : 2, 
                            "merged-rg-mkdup-iter-" + i + ".bam.gz" : 1, "merged-rg-mkdup-iter-" + i + ".bam.gz.bai" : 1, 
                            "pe-iter-" + i + ".bam.gz" : 2, "pem-iter-" + i + ".bam.gz" : 2, "se-iter-" + i + ".bam.gz" : 2 };

    possible_vcf_files = { "vcf-scaff" : 2, "gvcf-scaff" : 2, "iter-" + i + "-filter-intermediate.vcf.gz" : 2, "iter-" + i + "-filter-intermediate.vcf.gz.tbi" : 2, 
                            "iter-" + i + "-filter-intermediate-snps.vcf.gz" : 1, "iter-" + i + "-filter-intermediate-snps.vcf.gz.tbi" : 1, 
                            "iter-" + i + "-gathervcfs-params.txt" : 2, "iter-" + i + "-filter.vcf.gz" : 1, "iter-" + i + "-filter.vcf.gz.tbi" : 1, 
                            "iter-" + i + "-filter-snps.vcf.gz" : 1, "iter-" + i + "-filter-snps.vcf.gz.tbi" : 1 };

    # possible_fa_files = ["iter-" + prev_i + "-masked.fa", "iter-" + prev_i + "snps-masked.fa", "iter-" + i + "-snps-intermediate.dict", 
    #                         "iter-" + i + "-snps-intermediate.fa", "iter-" + i + "-snps-intermediate.fa.amb", "iter-" + i + "-snps-intermediate.fa.ann", 
    #                         "iter-" + i + "-snps-intermediate.fa.bwt", "iter-" + i + "-snps-intermediate.fa.fai", "iter-" + i + "-snps-intermediate.fa.pac", 
    #                         "iter-" + i + "-snps-intermediate.fa.sa" ];
    # All possible files/directories generated for each iteration. This also specifies their level in the keepall heirarchy: if a file's
    # level is above the current keep setting it is removed.
    # The list of fa files would need to be updated given the new suffixes above, but these are generally all left alone so
    # I just commented them out.

    ##########

    if globs['last-iter'] and globs['keeplevel'] == 0:
        globs['keeplevel'] = 1;
    # Keep more files for the last iteration.

    for f in possible_map_files:
        if possible_map_files[f] > globs['keeplevel']:
            full_f = os.path.join(globs['iterbamdir'], f);
            if os.path.isfile(full_f):
                cmd = "os.remove(" + full_f + ")";
                cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Removing file", 'outfile' : "", 'logfile' : "", 'start' : False };
                if globs['dryrun']:
                    PC.report_step(globs, cmds, cmd, "DRYRUN", cmd);
                else:
                    PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
                    os.remove(full_f);
    # Remove all mapping files above the keep threshold.

    ##########

    for f in possible_vcf_files:
        if possible_vcf_files[f] > globs['keeplevel']:
            full_f = os.path.join(globs['itervcfdir'], f);

            if os.path.isfile(full_f):
                cmd = "os.remove(" + full_f + ")";
                cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Removing file", 'outfile' : "", 'logfile' : "", 'start' : False };
                if globs['dryrun']:
                    PC.report_step(globs, cmds, cmd, "DRYRUN", cmd);
                else:    
                    PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
                    os.remove(full_f);
            elif os.path.isdir(full_f):
                cmd = "shutil.rmtree(" + full_f + ")";
                cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Removing directory", 'outfile' : "", 'logfile' : "", 'start' : False };
                if globs['dryrun']:
                    PC.report_step(globs, cmds, cmd, "DRYRUN", cmd);
                else:    
                    PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
                shutil.rmtree(full_f);
    # Remove all vcf files above the keep threshold.

    ##########

    return cmds;
    # Return the cmds dict with the remove commands we ran to clean up.

#############################################################################