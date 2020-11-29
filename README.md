# Pseudo-it
### Reference-guided genome assembly with iterative mapping

[![Python](https://goodest-goodlab.github.io/pseudo-it/img/pybadge.svg)](https://www.python.org/)
[![Build Status](https://travis-ci.org/goodest-goodlab/pseudo-it.svg?branch=master)](https://travis-ci.org/goodest-goodlab/pseudo-it)
[![codecov](https://codecov.io/gh/goodest-goodlab/pseudo-it/branch/master/graph/badge.svg?token=3U2K65R611)](https://codecov.io/gh/goodest-goodlab/pseudo-it)
[![Commit activity](https://img.shields.io/github/commit-activity/m/goodest-goodlab/pseudo-it)](https://github.com/goodest-goodlab/pseudo-it/commits/master)
[![License](https://img.shields.io/github/license/goodest-goodlab/pseudo-it)](https://github.com/goodest-goodlab/pseudo-it/LICENSE)

## Authors
#### Gregg Thomas, Brice Sarver, and Jeff Good

## About

### Pseudo-it is a program that uses an iterative read mapping process to generate pseudo-assemblies of closely related genomes to mitigate reference bias. 

#### What follows is a brief explanation of the options.

#### This is an update of version 2 of the software, found here: https://github.com/bricesarver/pseudo-it

## Citation

#### Sarver BAJ, Keeble S, Cosart T, Tucker PK, Dead MD, Good JM. 2017. Phylogenomic insights into mouse evolution using a pseudoreference approach. Genome Biology and Evololution. https://doi.org/10.1093/gbe/evx034.

## Version History
#### This is version beta 3.1, released August 13, 2020

Change log:
* Restructured code to add modularity for easy incorporation of other read mappers or variant callers.
* Now uses HaplotypeCaller from GATK4 for variant calling.
* Now uses bcftools to generate the final consensus FASTA file, which allows the incorporation of indels and generates a chain file to preserve coordinate system from original reference.

## Installation

#### Pseudo-it is implemented in Python 3.

Pseudo-it itself is simply a Python script that coordinates the running of other software.  Simply download the program and run it. You may want to add the Referee folder to your $PATH variable for ease of use.

Pseudo-it runs a suite of common bioinformatics software to go from raw reads in FASTQ format with a reference assembly in FASTA format to a pseudo-assembly of the reads in FASTA format. These programs are:

1. [BWA](http://bio-bwa.sourceforge.net/) for read mapping.
2. [GATK](https://gatk.broadinstitute.org/hc/en-us) for variant calling.
3. [samtools](http://www.htslib.org/download/) for handling mapped reads (BAM).
4. [Picard Tools](https://github.com/broadinstitute/picard) for handling mapped reads (BAM).
5. [bedtools](https://bedtools.readthedocs.io/en/latest/) for soft-masking FASTA files.
6. [bcftools](http://samtools.github.io/bcftools/) for handling VCF files and generating consensus FASTA and .chain files.

You will need each of these programs installed on your system to run pseudo-it. You can install them each individually, but they can also easily be installed from [bioconda](https://anaconda.org/bioconda), and we have provided the `pseudo-it.yml` pre-set environment to automatically install them.

First, make sure Anaconda is installed on your system by following the instructions here: [Anaconda Download](https://www.anaconda.com/products/individual#download-section)

Then, start Anaconda:
    
    source /path/to/anaconda3/bin/activate

Load the environment:

    conda env create --file psuedo-it.yml

And finally, activate the environment:

    conda activate pseudo-it

Additionally, Pseudo-it relies on common *nix based commands such as `grep`/`zgrep` and `sed`. **As such, Pseudo-it can only be run on systems with those commands installed**

### Verifying installation and dependencies

If you have added pseudo-it folder to your $PATH, you can run it from any location in your file system simply by typing `pseudo_it.py`.

Without adding it to your $PATH you will need to explicitly invoke python and provide the full path to the pseudo_it.py interface script:

    python /path/to/pseudo-it/pseudo_it.py

**All subsequent command examples assume you have added the pseudo-it folder to your $PATH**

#### Testing the pseudo-it interface:

To make sure pseudo_it can be called, simply use the version check:
    
    pseudo_it.py --version

This should print out something like this:

    # Pesudo-it version Beta 3.1 released on August 13, 2020

#### Checking that all dependencies are installed:

Option coming soon.

## Usage

The most basic usage of pseudo-it would be a command as follows:

    pseudo_it.py -ref [reference genome FASTA file] -pe1 [paired-end reads FASTQ file 1] -pe2 [paired-end reads FASTQ file 2] -i [number of iterations of mapping] -p [max number of processes to use] -o [desired output directory]

### Options

**Use `pseudo_it.py -h` to print out a help menu listing all the options.**

| Option | Description | 
| ------ | ----------- |
| -ref [FASTA file] | A FASTA formatted file containing the genome you wish to score. Can be gzip compressed or not. FASTA headers must match the sequence IDs in column one of the pileup or genotype log likelihood file. |
| -se [FASTQ file] | A FASTQ file containing single-end reads. |
| -pe1 [FASTQ file] | A FASTQ file containing pair 1 of paired-end reads. |
| -pe2 [FASTQ file] | A FASTQ file containing pair 2 of paired-end reads. | 
| -pem [FASTQ file] | A FASTQ file containing merged paired-end reads. |
| -bam [BAM file] | OPTIONAL: A BAM file with the reads provided in `-se`, `-pe1` and `-pe2`, and/or `-pem` already mapped to the reference assembly. |
| -f "[STRING]" | The expression to filter variants. Must conform to VCF INFO field standards. See [bcftools expressions](https://samtools.github.io/bcftools/bcftools.html#expressions) for more info. Default read depth filters are optimized for a 30-40X sequencing run -- adjust for your assembly. Default: "MQ < 30.0 || DP < 5 || DP > 60"
| -i [INT] | The number of iterations of pseudo-it to run. Default: 4 |
| -o [directory name] | The desired output directory. This will be created for you if it doesn't exist. Default: `pseudoit-out-[date]-[time]`. One of -o or -resume must be provided. |
| -resume [directory name] | The path to a previous Pseudo-it directory to resume a run. Scans for presence of files and resumes when it can't find an expected file. One of -o or -resume must be provided. |
| -tmp [directory name] | OPTIONAL: Some programs write files to a temporary directory. If your default tmp dir is size limited, specify a new one here, or just specifiy 'tmp-pi-out' to have a folder called 'tmp' created and used within the main output folder. |
| -p [INT] | The MAX number of processes Pseudo-it can use. It is highly recommended to run this with multiple processes. See [Resource allocation](#Resource-allocation) for more detail. |
| -bwa-t [INT] | The number of threads for BWA mem to use for each FASTQ library. See [Resource allocation](#Resource-allocation) for more detail. Default: 1. |
| -gatk-t [INT] | The number of threads for GATK's Haplotype caller to use. See [Resource allocation](#Resource-allocation) for more detail. Default: 4. |
| -bwa [PATH STRING] | The path to the BWA mapping progam. Default: bwa |
| -picard [PATH STRING] | The exact command used to run picard. For a jar file: `java -jar [full path to jar file]`. For an alias or conda install: picard. Include heap size in command if necessary, i.e. `-Xmx6g`. Default: picard |
| -samtools [PATH STRING] | The path to the samtools progam. Default: samtools |
| -gatk [PATH STRING] | The path to the GATK progam. Default: gatk |
| -bedtools [PATH STRING] | The path to the bedtools progam. Default: bedtools |
| -bcftools [PATH STRING] | The path to the bcftools progam. Default: bcftools |
| --maponly | Only do one iteration and stop after read mapping. |
| --noindels | Set this to not incorporate indels into the final assembly. |
| --diploid | Set this use IUPAC ambiguity codes in the final FASTA file. |
| --keepall | By default, pseudo-it keeps only the final files for each step of each iteration (BAM, VCF, FASTA and their respective indices). Set this option to keep all intermediate files. While this is the best way to ensure your runs can be resumed with different settings this will result many large files being saved (total of ~1TB for a 30X genome and 4 iterations). |
| --keeponlyfinal | By default, pseudo-it keeps only the final files for each step of each iteration (BAM, VCF, FASTA and their respective indices). Set this option to keep these files ONLY for the final iteration. While this minimizes storage space required, you will be unable to resume this run. |
| --overwrite | Set this to overwrite existing files (as opposed to -resume which skips steps that have files already writen). |
| --quiet | Set this flag to prevent psuedo-it from reporting detailed information about each step. |
| --version | Simply print the version and exit. Can also be called as `-version`, `-v`, or `--v`. |
| -h | Print a help menu and exit. Can also be called as `--help`. |

## Resource allocation

In addition to our own efforts to speed up the programs that pseudo-it runs, those programs themselves may also have the option to specify multiple processes/threads. This means that the allocation of processes/threads across programs can be confusing, and this section will try to explain how this is done.

1. First, multiple processes can be specified in the main pseudo-it interface with the `-p` option. This specifies the MAX number of processes to be used to run BWA and HaplotypeCaller. Minimally, each set of reads provided (`-se`, `-pe1` and `-pe2`, and/or `-pem`) would have 1 process for read mapping with BWA and each large scaffold/chromosome in the reference assembly would have 1 process for variant calling with HaplotypeCaller. 

2. BWA mem, the algorithm used to map reads by default, also allows for multiple threads with its internal `-t` option. Based on the number of read sets provided and the number of processes specified with `-p`, pseudo-it will automatically determine how many threads BWA mem should use.

    For example, if one set of paired end reads were supplied (`-pe1` and `-pe2`) and `-p` was set to 9, BWA mem would use all 9 threads for that single read set.

    If all three sets of reads are supplied (`-se`, `-pe1` and `-pe2`, and/or `-pem`) and `-p` was set to 9, BWA mem run all three sets simultaneously, each with 9 / 3 = 3 threads.

    Alternatively, you can specify the number of threads reserved for all BWA mem processes with the `-bwa-t` option. In the case of all three sets of reads being supplied with `-p 9` BUT with `-bwa-t 2`, each read set will be mapped simultaneously, but this time with only 2 threads a piece.

3. GATK's HaplotypeCaller also has a multi-thread option: `--native-pair-hmm-threads`. It is not clear how much this option effects runtime, but the default value with GATK is 4 which we leave as default for pseudo-it. 

    This means that if you run pseudo-it with a reference genome with 10 large scaffolds/chromosomes, you would ideally want to specify `-p 40` so that each chromosome gets its own GATK process that has access to 4 threads for the default `--native-pair-hmm-threads` option.

    Alternatively, with a reference genome with 10 large scaffolds/chromosomes, you could set pseudo-it's `gatk-t 1` to override the default `--native-pair-hmm-threads` option within GATK and you would only need `-p 10`.

    In the case of a fragmented assembly (i.e. 10 large scaffolds/chromosomes + 1000s of small scaffolds), you may want to specify one or two extra processes so the shorter scaffolds are not held up while the longer ones are being run.

### An optimized command for all 3 types of reads and a fragmented assembly with 10 large scaffolds + 1000s of small scaffolds

    python pseudo_it.py -ref [reference genome FASTA file] -se [single-end reads FASTQ file] -pe1 [paired-end reads FASTQ file 1] -pe2 [paired-end reads FASTQ file 2] -pem [merged paired-end reads FASTQ file] -i 4 -p 48 -o [desired output directory]

This command allocates 48 total processes for pseudo-it to use. It would run all three read sets through BWA mem with 48 / 3 = 16 threads each and would subsequently spawn 12 GATK HaplotypeCaller processes, each with the default `--native-pair-hmm-threads 4` option (12 * 4 = 48). The 10 long scaffolds would each have a separate processes while the shorter scaffolds would run have 2 extra processes allocated for them.

### An less resource-intensive command for all 3 types of reads and a fragmented assembly with 10 large scaffolds + 1000s of small scaffolds

    python pseudo_it.py -ref [reference genome FASTA file] -se [single-end reads FASTQ file] -pe1 [paired-end reads FASTQ file 1] -pe2 [paired-end reads FASTQ file 2] -pem [merged paired-end reads FASTQ file] -i 4 -p 12 -gatk-t 1 -o [desired output directory]

This command adds the `gatk-t 1` option to reduce the number of threads for HaplotypeCaller, thus reducing the total number of processes needed to `-p 12`. BWA mem would now run on each read set with only 4 threads.