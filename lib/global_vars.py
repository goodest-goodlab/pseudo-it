# This file holds some global variables for some of the input options.
# Global variables are exclusively read only -- they are not modified anywhere else in the code except when reading the input options.

import sys, timeit, lib.picore as PC

def init():
    globs = {
        'pyver' :  ".".join(map(str, sys.version_info[:3])),
        'version' : 'Beta 2.0',
        'releasedate' : 'December 07, 2019',
        'doi' : 'https://doi.org/10.1093/gbe/evx034',
        'http' : 'https://github.com/bricesarver/pseudo-it',
        'github' : 'https://github.com/bricesarver/pseudo-it/issues',
        'starttime' : timeit.default_timer(),
        'startdatetime' : PC.getOutTime(),
        # Meta info

        'ref' : False,
        'scaffs' : False,
        # Reference FASTA file

        'se' : False,
        'pe1' : False,
        'pe2' : False,
        'pem' : False,
        'libs' : {},
        'num-libs' : 0,
        # FASTQ library info

        'indir' : '',
        'outdir' : '',
        'logfilename' : 'pseudo-it.errlog',
        'logdir' : '',
        'tmpdir' : 'System default.',
        # I/O stuff
        
        'bwa-path' : 'bwa',
        'picard-path' : 'picard',
        'samtools-path' : 'samtools',
        'gatk-path' : 'gatk',
        'bedtools-path' : 'bedtools',
        'bcftools-path' : 'bcftools',
        # Dependency paths

        'bwa-t' : 1,
        # Number of threads for BWA mem to use.

        'gatk-t' : 4,
        # Number of threads for GATK's --native-pair-hmm-threads option.     

        'heap' : 'System default.',
        # Heap size for java programs (picard).

        'filter' : '"MQ < 30.0 || FORMAT/DP < 5 || FORMAT/DP > 60"',
        # Variant filtration string default

        'indels' : True,
        'resume' : False,
        'overwrite' : False,
        'num-procs' : 1,
        'num-iters' : 4,
        'quiet' : False,
        # User options

        'norun' : False,
        'dryrun' : False,
        'continue' : False,
        'iteration' : 1,
        'pad' : 82,
        'endprog' : False,
        'log-v' : 1,
        'stats' : True,
        'progstarttime' : 0,
        'stepstarttime' : 0,
        'pids' : "",
        'psutil' : "",
        'debug' : False,
        'nolog' : False,
        # Internal stuff
    }

    globs['logfilename'] = "pseudo-it-" + globs['startdatetime'] + ".errlog";
    # Add the runtime to the error log file name.

    return globs;