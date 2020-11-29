# Consensus reference generation for Pseudo-it
#############################################################################
import os, lib.picore as PC, lib.piref as piref
#############################################################################

def getMask(globs, cmds, vcf_file):
# Get the sites to be masked into a bed file.

    mask_bedfile = os.path.join(globs['iterfadir'], "iter-" + globs['iter-str'] + "-masksites.bed");
    if globs['diploid']:
        mask_bedfile = mask_bedfile.replace("-masksites.bed", "-diploid-masksites.bed");

    cmd = "zgrep \"\./\.\" " + vcf_file + " | awk '{{OFS=\"\t\"; if ($0 !~ /\#/); print $1, $2-1, $2}}' | bedtools merge -i - > " + mask_bedfile;
    cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Get mask sites", 'outfile' : mask_bedfile, 'logfile' : "", 'start' : False };

    run = True;
    if globs['resume']:
        if os.path.isfile(mask_bedfile) and os.stat(mask_bedfile).st_size != 0:
            PC.report_step(globs, cmds, cmd, "RESUME", "previous output found: " + mask_bedfile);
            run = False;
            
    if run:
        if not globs['dryrun']:
            PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
            os.system(cmd);

            if os.path.isfile(mask_bedfile) and os.stat(mask_bedfile).st_size != 0:
                num_sites = str(len(open(mask_bedfile, "r").readlines()));
                PC.report_step(globs, cmds, cmd, "SUCCESS", num_sites + " mask sites read: " + mask_bedfile);
            else:
                PC.report_step(globs, cmds, cmd, "ERROR!", "Mask sites file not found or empty: " + mask_bedfile);
                globs['exit-code'] = 1;
                PC.endProg(globs);

        else:
            PC.report_step(globs, cmds, cmd, "DRYRUN", cmd);

    return mask_bedfile, cmds; 

#############################################################################

def maskFa(globs, cmds, mask_bedfile, cur_ref):
# Fix the headers from the consensus FASTA file.
    prev_iter = str(int(globs['iter-str']) - 1);
    if len(prev_iter) == 1:
        prev_iter = "0" + prev_iter;

    if globs['indels']:
        cur_logfile = os.path.join(globs['iterlogdir'], "bedtools-maskfasta-" + globs['iter-str'] + ".log");
        mask_ref = os.path.join(globs['iterfadir'], "iter-" + prev_iter + "-masked.fa");
    else:
        cur_logfile = os.path.join(globs['iterlogdir'], "bedtools-maskfasta-" + globs['iter-str'] + "-snps.log");
        mask_ref = os.path.join(globs['iterfadir'], "iter-" + prev_iter + "-snps-masked.fa");

    mask_cmd = globs['bedtools-path'] + " maskfasta -fi " + cur_ref + " -bed " + mask_bedfile + " -soft -fo " + mask_ref;
    cmds[mask_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Softmask reference", 'outfile' : mask_ref, 'logfile' : cur_logfile, 'start' : False };
    
    exit_flag = PC.runCMD(mask_cmd, globs, cmds, True);
    PC.exitCheck(exit_flag, globs);
    # End the program if an error is encountered  

    return mask_ref, cmds;

#############################################################################

def genConsensus(globs, cmds, vcf_file, cur_ref):
# Run the command to generate a consensus FASTA file from the reference and the variants.
    
    cmd = "getConsCase()";
    cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Determining case of first base", 'outfile' : "", 'logfile' : "", 'start' : False };
    
    bcftools_cmd = globs['bcftools-path'] + " consensus -f " + cur_ref + " -o " + globs['iter-final-fa']
    if globs['last-iter'] and globs['indels']:
        bcftools_cmd += " -c " + globs['iter-final-chain'];
    if globs['last-iter'] and globs['diploid']:
        bcftools_cmd += " -I ";
    bcftools_cmd += " -e \"FILTER='pseudoit'\" " + vcf_file;
    cmds[bcftools_cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Generating consensus", 'outfile' : globs['iter-final-fa'], 'logfile' : globs['iter-consensus-log'], 'start' : False };

    run_flag = True;
    if globs['resume']:
        run_flag = PC.runCheck(bcftools_cmd, cmds, globs);

    #### RUN RUNCHECK FIRST

    first_lower = False;
    if globs['dryrun']:
        PC.report_step(globs, cmds, cmd, "DRYRUN", cmd);
        first_lower, linestr_orig, linestr_repl = True, "a", "A";
    elif run_flag:   
        PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
        first_lower, linestr_orig, linestr_repl = getConsCase(cur_ref);
        PC.report_step(globs, cmds, cmd, "SUCCESS", "First base: " + linestr_orig[0]);
    # This first_lower stuff is a hack to deal with bcftools consensus using the case of the first base in the reference fasta to inject variants.
    # Possibly resolved: https://github.com/samtools/bcftools/issues/1150#issuecomment-582407490
    # Need to test and make sure it is in official release before I remove this hack.

    if first_lower:
        cmd = "sed -i '2 s/" + linestr_orig + "/" + linestr_repl + "/g' " + cur_ref;
        cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Changing first ref base to upper case", 'outfile' : "", 'logfile' : "", 'start' : False };

        if globs['dryrun']:
            PC.report_step(globs, cmds, cmd, "DRYRUN", cmd);
        elif run_flag:  
            PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
            os.system(cmd);
            PC.report_step(globs, cmds, cmd, "SUCCESS", "First base converted to upper case");
    # Part of first_lower hack.

    exit_flag = PC.runCMD(bcftools_cmd, globs, cmds, True);    
    # Consensus command
    PC.exitCheck(exit_flag, globs);
    # End the program if an error is encountered  

    if first_lower:
        cmd = "sed -i '2 s/" + linestr_repl + "/" + linestr_orig + "/g' " + cur_ref;
        cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Reverting case of first ref base", 'outfile' : "", 'logfile' : "", 'start' : False };
        
        if globs['dryrun']:
            PC.report_step(globs, cmds, cmd, "DRYRUN", cmd);
        elif run_flag: 
            PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
            os.system(cmd);
            PC.report_step(globs, cmds, cmd, "SUCCESS", "First base reverted to original case");


        if not globs['dryrun']:
            first_lower, linestr_orig, linestr_repl = getConsCase(globs['iter-final-fa']);

        cmd = "sed -i '2 s/" + linestr_orig + "/" + linestr_repl + "/g' " + globs['iter-final-fa'];
        cmds[cmd] = { 'cmd-num' : PC.getCMDNum(globs, len(cmds)), 'desc' : "Reverting case of first consensus base", 'outfile' : "", 'logfile' : "", 'start' : False };
        
        if globs['dryrun']:
            PC.report_step(globs, cmds, cmd, "DRYRUN", cmd);
        elif run_flag:            
            PC.report_step(globs, cmds, cmd, "EXECUTING", cmd);
            os.system(cmd);
            PC.report_step(globs, cmds, cmd, "SUCCESS", "First base reverted to original case");
    # Part of first_lower hack.

    globs['consensus-file'] = globs['iter-final-fa'];

    return cmds, globs;

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