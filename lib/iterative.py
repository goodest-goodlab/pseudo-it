# Iterative mapping steps for Pseudo-it
#############################################################################
import sys, os, shutil, subprocess, lib.picore as PC, lib.piref as piref, \
    lib.pimap as pimap, multiprocessing as mp, \
    lib.varcall as varcall, lib.varpost as varpost, lib.consensus as con
#############################################################################

def mapping(globs):

    globs = PC.getIterStr(globs);
    PC.printWrite(globs['logfilename'], globs['log-v'], "#\n# " + "=" * 51 + " ITERATION " + globs['iter-str'] + " STARTING! " + "=" * 50);

    globs['iterstarttime'] = PC.report_step(globs, "", start=True);
    if globs['iteration'] == 1:
        globs['progstarttime'] = globs['iterstarttime'];

    cmds = {};
    
    cur_ref = PC.getRef(globs);
    globs['last-iter'] = False;
    if globs['iteration'] == globs['num-iters']:
        globs['last-iter'] = True;
    # Iteration prep

    globs['iterdir'] = os.path.join(globs['outdir'], "iter-" + globs['iter-str']);
    globs['iterbamdir'] = os.path.join(globs['iterdir'], "bam");
    globs['itervcfdir'] = os.path.join(globs['iterdir'], "vcf");
    globs['itervcfdir'] = os.path.join(globs['iterdir'], "vcf");
    globs['iterfadir'] = os.path.join(globs['iterdir'], "fa");
    globs['iterlogdir'] = os.path.join(globs['iterdir'], "logs");
    for d in [ globs['iterdir'], globs['iterbamdir'], globs['itervcfdir'], globs['iterfadir'], globs['iterlogdir'] ]:
        if not os.path.isdir(d):
            os.makedirs(d);
    if globs['last-iter']:
        globs['itervcfscaffdir'] = os.path.join(globs['itervcfdir'], "gvcf-scaff");
        globs['itervcflogdir'] = os.path.join(globs['itervcfdir'], "gvcf-logs");
    else:
        globs['itervcfscaffdir'] = os.path.join(globs['itervcfdir'], "vcf-scaff");
        globs['itervcflogdir'] = os.path.join(globs['itervcfdir'], "vcf-logs");
    for d in [ globs['itervcfscaffdir'], globs['itervcflogdir'] ]:
        if not os.path.isdir(d):
            os.makedirs(d);
    # Make directories for current iteration

    if globs['bam'] and globs['iteration'] == 1:
        globs['iter-final-bam-log'] = "NA";
        globs['iter-final-bam'] = globs['bam'];  
    else:      
        globs['iter-final-bam-log'] = os.path.join(globs['iterlogdir'], "picard-mkdup-iter-" + globs['iter-str'] + ".log");
        globs['iter-final-bam'] = os.path.join(globs['iterbamdir'], "merged-rg-mkdup-iter-" + globs['iter-str'] + ".bam.gz");
    # Final BAM file for this iteration

    if globs['last-iter']:
        globs['iter-gather-vcf-log'] = os.path.join(globs['iterlogdir'], "gatk-gathervcfs-iter-" + globs['iter-str'] + ".log");
        globs['iter-gather-vcf'] = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + "-filter.vcf.gz");

        if globs['indels']:
            globs['iter-final-vcf-log'] = globs['iter-gather-vcf-log'];
            globs['iter-final-vcf'] = globs['iter-gather-vcf'];
        else:
            globs['iter-final-vcf-log'] = os.path.join(globs['iterlogdir'], "gatk-selectsnps-iter-" + globs['iter-str'] + ".log");
            globs['iter-final-vcf'] = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + "-filter-snps.vcf.gz");
    else:
        globs['iter-gather-vcf-log'] = os.path.join(globs['iterlogdir'], "gatk-gathervcfs-iter-" + globs['iter-str'] + "-intermediate.log");
        globs['iter-gather-vcf'] = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + "-filter-intermediate.vcf.gz");

        globs['iter-final-vcf-log'] = os.path.join(globs['iterlogdir'], "gatk-selectsnps-iter-" + globs['iter-str'] + "-intermediate.log");
        globs['iter-final-vcf'] = os.path.join(globs['itervcfdir'], "iter-" + globs['iter-str'] + "-filter-intermediate-snps.vcf.gz");
    # Final VCF file for this iteration

    if globs['last-iter']:
        if globs['indels']:
            globs['iter-consensus-log'] = os.path.join(globs['iterlogdir'], "bcftools-consensus-iter-" + globs['iter-str'] + "-final.log");
            globs['iter-final-chain'] = os.path.join(globs['iterfadir'], "iter-" + globs['iter-str'] + "-final.chain");
            globs['iter-final-fa'] = os.path.join(globs['iterfadir'], "iter-" + globs['iter-str'] + "-final.fa");

        else:
            globs['iter-consensus-log'] = os.path.join(globs['iterlogdir'], "bcftools-consensus-iter-" + globs['iter-str'] + "-snps-final.log");
            globs['iter-final-chain'] = os.path.join(globs['iterfadir'], "iter-" + globs['iter-str'] + "-snps-final.chain");
            globs['iter-final-fa'] = os.path.join(globs['iterfadir'], "iter-" + globs['iter-str'] + "-snps-final.fa");

        if globs['diploid']:
            globs['iter-consensus-log'] = globs['iter-consensus-log'].replace("-final.log", "-diploid-final.log");
            globs['iter-final-chain'] = globs['iter-final-chain'].replace("-final.chain", "-diploid-final.chain");
            globs['iter-final-fa'] = globs['iter-final-fa'].replace("-final.fa", "-diploid-final.fa");

    else:
        globs['iter-consensus-log'] = os.path.join(globs['iterlogdir'], "bcftools-consensus-iter-" + globs['iter-str'] + "-snps-intermediate.log");
        globs['iter-final-chain'] = os.path.join(globs['iterfadir'], "iter-" + globs['iter-str'] + "-snps-intermediate.chain");
        globs['iter-final-fa'] = os.path.join(globs['iterfadir'], "iter-" + globs['iter-str'] + "-snps-intermediate.fa");
    # Final FASTA files for this iteration
    # Output files for current iteration

    if globs['iteration'] == 1:
        cmds = piref.indexCheck(globs['ref'], globs, cmds);
        cmds = piref.getScaffs(globs['ref'], globs, cmds);
    # Check that all index files have been created and get the scaffold IDs from the reference FASTA

    if globs['bam'] and globs['iteration'] == 1:
        do_mapping = False;
    elif globs['resume']:
        do_mapping = PC.prevCheck(globs['iter-final-bam'], globs['iter-final-bam-log'], globs);
    else:
        do_mapping = True;

    if do_mapping:
        do_varcalling = True;
    elif globs['resume']:
        do_varcalling = PC.prevCheck(globs['iter-final-vcf'], globs['iter-final-vcf-log'], globs);
    else:
        do_varcalling = True;

    if do_varcalling:
        do_consensus = True;
    elif globs['resume']:
        do_consensus = PC.prevCheck(globs['iter-final-fa'], globs['iter-consensus-log'], globs);
    else:
        do_consensus = True;
    # CHECK WHICH STEPS WE NEED TO PERFORM

    statstr = "EXECUTING";
    if globs['resume']:
        statstr = "RESUME";
    if globs['dryrun']:
        statstr = "DRYRUN";
    # Status for the main step reports

    if globs['iteration'] != 1:
        cmds = piref.indexFa(globs, cmds, cur_ref);
    # INDEX FASTA IF NOT FIRST ITERATION

    if do_mapping:
        PC.report_step(globs, cmds, "NA--01   Read mapping", statstr, "Mapping reads and post-processing.");
        bamfiles, cmds = pimap.BWA(globs, cmds, cur_ref);
        # READ MAPPING

        #rg_bamfile, cmds = varprep.addRG(globs, cmds, bamfiles);
        # ADD READ GROUPS

        merged_bamfile, cmds = pimap.mergeBam(globs, cmds, bamfiles);
        # MERGE BAM FILES also sorts

        cmds = pimap.markDups(globs, cmds, merged_bamfile);
        # MARK DUPLICATES
    elif globs['bam'] and globs['iteration'] == 1:
        PC.report_step(globs, cmds, "NA--01   Read mapping", "BAM", "initial BAM file provided, skipping all mapping steps: " + globs['iter-final-bam']);
    else:
        PC.report_step(globs, cmds, "NA--01   Read mapping", "RESUME", "previous processed BAM file found, skipping all mapping steps: " + globs['iter-final-bam']);
    ## READ MAPPING STEPS

    #cmds = varprep.indexBAM(globs, cmds);
    # INDEX BAM FILE
    # Now done during MarkDuplicates

    if globs['map-only']:
        PC.report_step(globs, cmds, "NA--04   Iteration cleanup", statstr, "Removing intermediate files based on --keep* options.");
        cmds = cleanUp(globs, cmds);

        PC.printWrite(globs['logfilename'], globs['log-v'], "#\n# " + "=" * 51 + " ITERATION " + globs['iter-str'] + " COMPLETE! " + "=" * 50);
        PC.report_step(globs, "", end=True);

        globs['iteration'] +=1;
        return globs;
    # This stops the program after the first iteration of mapping if --maponly is set.


    if do_varcalling:
        PC.report_step(globs, cmds, "NA--02   Variant calling", statstr, "Calling and post-processing variants.");
        cmds = varcall.haplotypeCaller(globs, cmds, cur_ref, globs['iter-final-bam']);
        # HAPLOTYPECALLER

        if globs['last-iter']:
            cmds = varcall.genotypeGVCFs(globs, cmds, cur_ref);
        # GENOTYPE GVCFS FOR LAST ITER

        # if not globs['last-iter'] or (globs['last-iter'] and not globs['indels']):
        #     cmds = varpost.selectSNPs(globs, cmds);
        # SELECT SNPs if it is not the last iteration, or if it is and the final output should not contain indels

        cmds = varpost.varFilter(globs, cmds, cur_ref);
        # FILTER VCFs

        cmds = varpost.gatherVCFs(globs, cmds);
        # COMBINE VCF

        cmds = varpost.indexVCF(globs, cmds, globs['iter-gather-vcf']);
        # INDEX VCF        
    else:
        PC.report_step(globs, cmds, "NA--02   Variant calling", "RESUME", "previous processed VCF file found, skipping all variant calling steps: " + globs['iter-final-vcf']);
    ## VARIANT CALLING STEPS
   
    if do_consensus:
        PC.report_step(globs, cmds, "NA--03   Consensus generation", statstr, "Generating consensus FASTA.");
        if globs['last-iter']:
            mask_bedfile, cmds = con.getMask(globs, cmds, globs['iter-gather-vcf']);
            # GET MASK SITES

            cur_ref, cmds = con.maskFa(globs, cmds, mask_bedfile, cur_ref);
            # MASK PREVIOUS REFERENCE

        if not globs['last-iter'] or (globs['last-iter'] and not globs['indels']):
            cmds = varpost.selectSNPs(globs, cmds, globs['iter-gather-vcf']);
            # SELECT SNPs FROM VCF IF IT IS NOT THE LAST ITERATION OR IF --noindels IS SET

            cmds = varpost.indexVCF(globs, cmds, globs['iter-final-vcf']);
            # INDEX FINAL VCF

        cmds, globs = con.genConsensus(globs, cmds, globs['iter-final-vcf'], cur_ref);
        # GENERATE CONSENSUS
    else:
        PC.report_step(globs, cmds, "NA--03   Consensus generation", "RESUME", "previous processed consensus FASTA file found, skipping all consensus generation steps: " + globs['iter-final-fa']);    
        globs['consensus-file'] = globs['iter-final-fa'];
    ## CONSENSUS STEPS


    PC.report_step(globs, cmds, "NA--04   Iteration cleanup", statstr, "Removing intermediate files based on --keep* options.");
    cmds = cleanUp(globs, cmds);

    PC.printWrite(globs['logfilename'], globs['log-v'], "#\n# " + "=" * 51 + " ITERATION " + globs['iter-str'] + " COMPLETE! " + "=" * 50);
    PC.report_step(globs, "", end=True);

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

    possible_fa_files = ["iter-" + prev_i + "-masked.fa", "iter-" + prev_i + "snps-masked.fa", "iter-" + i + "-snps-intermediate.dict", 
                            "iter-" + i + "-snps-intermediate.fa", "iter-" + i + "-snps-intermediate.fa.amb", "iter-" + i + "-snps-intermediate.fa.ann", 
                            "iter-" + i + "-snps-intermediate.fa.bwt", "iter-" + i + "-snps-intermediate.fa.fai", "iter-" + i + "-snps-intermediate.fa.pac", 
                            "iter-" + i + "-snps-intermediate.fa.sa" ]


    if globs['last-iter'] and globs['keeplevel'] == 0:
        globs['keeplevel'] = 1;

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


    return cmds;