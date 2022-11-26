# This file holds some global variables for some of the input options.
# Global variables are exclusively read only -- they are not modified anywhere else in the code except when reading the input options.

import sys
import timeit
import pseudo_it_lib.picore as PC

def init():
    globs = {
        
        'version' : 'Beta 3.1',
        'releasedate' : 'August 13, 2020',
        'doi' : 'https://doi.org/10.1093/gbe/evx034',
        'http' : 'https://github.com/bricesarver/pseudo-it',
        'github' : 'https://github.com/bricesarver/pseudo-it/issues',
        'starttime' : timeit.default_timer(),
        'startdatetime' : PC.getOutTime(),
        # Meta info

        'pyver' :  ".".join(map(str, sys.version_info[:3])),
        # System info

        'call' : "",
        # Script call info

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

        'bam' : False,
        'bam-index' : False,
        'indir' : '',
        'outdir' : '',
        'sample-name' : '',
        'logfilename' : 'pseudo-it.errlog',
        'logdir' : '',
        'tmpdir' : 'System default.',
        'indels' : True,
        'diploid' : False,
        'resume' : False,
        'overwrite' : False,
        'final' : False,
        'keeplevel' : 1,
        'mkdups' : True,
        'mask' : 'soft',
        'in-vcf' : False,
        'filter-sites' : False,
        # I/O options
        
        'mapper' : 'bwa',
        'mapper-path' : 'bwa',
        'picard-path' : 'picard',
        'samtools-path' : 'samtools',
        'gatk-path' : 'gatk',
        'bedtools-path' : 'bedtools',
        'bcftools-path' : 'bcftools',
        'tabix-path' : 'tabix',
        # Dependency paths

        'mapper-t' : 1,
        # Number of threads for the mapper to use.

        'gatk-t' : 4,
        'gvcf-procs' : 1,
        'filter-procs' : 1,
        # Number of threads for GATK's --native-pair-hmm-threads option and the number of procs
        # to use for GenotypeGVCFs and bcftools filter.

        'rg' : { 'ID' : '',
                 'PL' : '',
                 'PU' : '',
                 'LB' : '',
                 'SM' : ''},
        # Read groups.

        'filter' : "'MQ < 30.0 || FORMAT/DP < 5 || FORMAT/DP > 60'",
        # Variant filtration string default

        'in-bed' : False,
        'regions' : False,
        'bed-mode' : False,
        # Variables if a bed file containing intervals is provided. Variants will only be called in these intervals.
        # bed-mode can be one of 'regions' or 'file'.

        'intervals' : False,
        # The intervales over which to parallelize (or not) variant calling. 
        # By default, the scaffolds/chromosomes in the reference genome. 
        # If a bed file is specified with -bedmode regions, then those regions in the bed file. This is slow for many small regions but better for large regions.
        # If a bed file is specifed with -bedmode file, then a single gatk instance will spawn for the whole file with all specified processes. Better for many small regions.

        'num-procs' : 1,
        'num-iters' : 4,
        'quiet' : False,
        # Other user options

        'norun' : False,
        'dryrun' : False,
        'continue' : False,
        'iteration' : 1,
        'last-iter' : False,
        'map-only' : False,
        'pad' : 82,
        'endprog' : False,
        'exit-code' : 0,
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