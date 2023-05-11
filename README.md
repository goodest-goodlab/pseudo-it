# <p align="center"><img align="center" width="360" height="115" src="https://goodest-goodlab.github.io/pseudo-it/img/logo01.png"></p>

## Pseudo-it: Reference-guided genome assembly with iterative mapping

[![OS](https://goodest-goodlab.github.io/pseudo-it/img/osbadge.svg)](#Installation)
[![Language](https://goodest-goodlab.github.io/pseudo-it/img/pybadge.svg)](https://www.python.org/)
[![Build Status](https://travis-ci.org/goodest-goodlab/pseudo-it.svg?branch=master)](https://travis-ci.org/goodest-goodlab/pseudo-it)
[![codecov](https://codecov.io/gh/goodest-goodlab/pseudo-it/branch/master/graph/badge.svg?token=3U2K65R611)](https://codecov.io/gh/goodest-goodlab/pseudo-it)
![GitHub last commit](https://img.shields.io/github/last-commit/goodest-goodlab/pseudo-it)
[![License](https://img.shields.io/github/license/goodest-goodlab/pseudo-it)](https://github.com/goodest-goodlab/pseudo-it/LICENSE)

# Authors
### Gregg Thomas, Brice Sarver, and Jeff Good

# Table of Contents

- [About](#about)
    - [Citation](#citation)
- [Installation](#installation)
    - [Setting up conda](#0-install-and-set-up-conda-required-for-both-options-below)
    - [Installing pseudo-it from bioconda](#1-install-pseudo-it-through-bioconda-recommended)
    - [Installing pseudo-it from github and dependencies from conda](#2-install-pseudo-it-from-github-and-dependencies-with-conda)
    - [Checking that all dependencies are installed](#checking-that-all-dependencies-are-installed)
- [Usage](#usage)
    - [Options](#options)
    - [Pre-indexing your reference genome](#pre-indexing-your-reference-genome)
    - [Resource allocation](#resource-allocation)
    - [Example commands](#example-commands)
- [FAQs](#faqs)

# About

### Pseudo-it is a program that uses an iterative read mapping process to generate pseudo-assemblies of closely related genomes to mitigate reference bias. 

Pseudo-it performs a set of read mapping and variant calling steps to go from raw reads in FASTQ format with reference genome in FASTA format to a pseudo-assembly of the provided reads in FASTA format. This process is done iteratively to increase the number of reads mapped and capture more variation in the final assembly.

For the first iteration, pseudo-it performs the following steps:

1. Map provided reads (FASTQ) to **provided reference genome (FASTA).**
2. Call variants on mapped reads.
3. Generate a consensus FASTA file by inserting the called variants into the original sequence.

For each subsequent iteration, the previous iteration's consensus FASTA file serves as the new reference for read mapping:

1. Map provided reads (FASTQ) to **previous iteration consensus sequence (FASTA).**
2. Call variants on mapped reads.
3. Generate a new consensus FASTA file by inserting the called variants into the previous iteration's sequence.

Each iteration should allow for more reads to be mapped and more variation to be incorporated into the assembly.

#### What follows is a brief explanation of the options.

#### This is an update of version 2 of the software, found here: https://github.com/bricesarver/pseudo-it

## Citation

#### Sarver BAJ, Keeble S, Cosart T, Tucker PK, Dead MD, Good JM. 2017. Phylogenomic insights into mouse evolution using a pseudoreference approach. Genome Biology and Evololution. https://doi.org/10.1093/gbe/evx034.

# Installation

#### Pseudo-it is implemented in Python 3.

Pseudo-it itself is simply a Python script that coordinates the running of other software. As such, the only dependency of the script is Python 3, however pseudo-it runs a suite of common bioinformatics software to go from raw reads in FASTQ format with a reference assembly in FASTA format to a pseudo-assembly of the reads in FASTA format. These programs are:

1. [BWA](http://bio-bwa.sourceforge.net/) for read mapping.
2. [GATK](https://gatk.broadinstitute.org/hc/en-us) for variant calling.
3. [samtools](http://www.htslib.org/download/) for handling mapped reads (BAM).
4. [Picard Tools](https://github.com/broadinstitute/picard) for handling mapped reads (BAM).
5. [bedtools](https://bedtools.readthedocs.io/en/latest/) for soft-masking FASTA files.
6. [bcftools](http://samtools.github.io/bcftools/) for handling VCF files and generating consensus FASTA and .chain files.

Additionally, Pseudo-it relies on common *nix based commands such as `grep`/`zgrep` and `sed`. **As such, Pseudo-it can only be run on systems with those commands installed**

## Two ways to install pseudo-it and its dependencies

There are two ways to install pseudo-it and its dependencies. BOTH require the use of the [conda package manager](https://docs.conda.io/en/latest/) or the speedier solver, [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html)

### 0. Install and set-up conda (REQUIRED for both options below)

To set these up, do the following:

1. (Install conda)[https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html]

2. (Optional) (Install mamba)[https://mamba.readthedocs.io/en/latest/installation.html] -- note that if you do this step, all subsequent commands where `conda` is used can be replaced with `mamba`.

3. (Set-up your conda channels for bioconda)[https://bioconda.github.io/#usage]

### 1. Install pseudo-it through bioconda (RECOMMENDED)

With conda activated,

1. Create a new enviornment for pseudo-it:

```bash
conda create -n pseudo-it-env
```

2. Activate the `pseudo-it-env` environment:

```bash
conda activate pseudo-it-env
```

3. Install the pseudo-it package:


```bash
conda install pseudo-it
```

This will install the pseudo-it script as well as all the dependencies listed above. Now you can run pseudo-it as follows:

```bash
python pseudo_it.py -h
```

### 2. Install pseudo-it from github and dependencies with conda

1. Since pseudo-it itself is just a python script, you can download it directly from this github either by clicking the green **Code** button above and copying the link to `git clone` the repo, or by downloading the release on the right. Once downloaded, you may want to add this folder to your `$PATH` variable so you don't have to type the full path to the script every time you run it.

2. For the dependencies, we have provided an environment file (`pseudo-it.yml`) that will automatically download and install them from bioconda. First, load the environment:

```bash
conda env create --file pseudo-it.yml
```

3. Activate the newly created environment:

```bash
conda activate pseudo-it
```

4. Confirm the dependencies have been installed:

```bash
python /PATH/TO/pseudo_it.py --depcheck
```

Here, the path to the pseudo_it.py script is wherever you downloaded it to in step 1. If you added the pseudo-it folder to your `$PATH` variable, you can just call it as `python pseudo_it.py`.

**All subsequent command examples assume you have added the pseudo-it folder to your $PATH**

## Verifying installation and dependencies

### Testing the pseudo-it interface:

To make sure pseudo_it can be called, simply use the version check:
    
    pseudo_it.py --version

This should print out something like this:

    # Pesudo-it version Beta 3.1 released on August 13, 2020

### Checking that all dependencies are installed:

To make sure all dependencies have been properly installed and the correct paths given, use the `--depcheck` option:

    pseudo_it.py --depcheck

This should try to run each external program from within pseudo-it and print whether it was able to successfully execute each one. If everything works, you should see a nice table:

    # --depcheck set: CHECKING DEPENDENCY PATHS AND EXITING.

    PROGRAM    PATH          STATUS
    -------------------------------
    BWA        bwa           PASSED
    Picard     picard        PASSED
    samtools   samtools      PASSED
    GATK       gatk          PASSED
    bedtools   bedtools      PASSED
    bcftools   bcftools      PASSED
    tabix      tabix         PASSED

    # All dependencies PASSED.

# Usage

The most basic usage of pseudo-it would be a command as follows:

    pseudo_it.py -ref [reference genome FASTA file] -pe1 [paired-end reads FASTQ file 1] -pe2 [paired-end reads FASTQ file 2] -i [number of iterations of mapping] -p [max number of processes to use] -o [desired output directory]

## Options

**Use `pseudo_it.py -h` to print out a help menu listing all the options.**

Note that only one of `-se`, `-pe1` and `-pe2`, or `-pem` is required, but any combination of the three read types is acceptable.

| Option | Description | 
| ------ | ----------- |
| `-ref [FASTA file]` | A FASTA formatted file containing the genome you wish to map reads to in the first iteration. |
| `-se [FASTQ file]` | A FASTQ file containing single-end reads. |
| `-pe1 [FASTQ file]` | A FASTQ file containing pair 1 of paired-end reads. |
| `-pe2 [FASTQ file]` | A FASTQ file containing pair 2 of paired-end reads. | 
| `-pem [FASTQ file]` | A FASTQ file containing merged paired-end reads. |
| `-bam [BAM file]` | OPTIONAL: A BAM file with the reads provided in `-se`, `-pe1` and `-pe2`, and/or `-pem` already mapped to the reference assembly. |
| `-f "[STRING]"` | The expression to filter variants. Must conform to VCF INFO field standards. See [bcftools expressions](https://samtools.github.io/bcftools/bcftools.html#expressions) for more info. Default read depth filters are optimized for a 30-40X sequencing run -- adjust for your assembly. Default: `"MQ < 30.0 \|\| DP < 5 \|\| DP > 60"`
| `-i [INT]` | The number of iterations of pseudo-it to run. Default: 4 |
| `-o [directory name]` | The desired output directory. This will be created for you if it doesn't exist. Default: `pseudoit-out-[date]-[time]`. One of `-o` or `-resume` must be provided. |
| `-resume [directory name]` | The path to a previous Pseudo-it directory to resume a run. Scans for presence of files and resumes when it can't find an expected file. One of `-o` or `-resume` must be provided. |
| `-tmp [directory name]` | OPTIONAL: Some programs write files to a temporary directory. If your default tmp dir is size limited, specify a new one here, or just specifiy 'tmp-pi-out' to have a folder called 'tmp' created and used within the main output folder. |
| `-p [INT]` | The MAX number of processes Pseudo-it can use. It is highly recommended to run this with multiple processes. See [Resource allocation](#Resource-allocation) for more detail. |
| `-bwa-t [INT]` | The number of threads for BWA mem to use for each FASTQ library. See [Resource allocation](#Resource-allocation) for more detail. Default: 1. |
| `-gatk-t [INT]` | The number of threads for GATK's Haplotype caller to use. See [Resource allocation](#Resource-allocation) for more detail. Default: 4. |
| `-bwa [PATH STRING]` | The path to the BWA mapping progam. Default: bwa |
| `-picard [PATH STRING]` | The exact command used to run picard. For a jar file: `java -jar [full path to jar file]`. For an alias or conda install: picard. Include heap size in command if necessary, i.e. `-Xmx6g`. Default: picard |
| `-samtools [PATH STRING]` | The path to the samtools progam. Default: samtools |
| `-gatk [PATH STRING]` | The path to the GATK progam. Default: gatk |
| `-bedtools [PATH STRING]` | The path to the bedtools progam. Default: bedtools |
| `-bcftools [PATH STRING]` | The path to the bcftools progam. Default: bcftools |
| `--depcheck` | Check that all dependencies can be executed from within pseudo-it. Paths can be provided with the option for each program, or left as default. |
| `--dryrun` | Run through all of the commands that would be run with the given inputs. Good to run to setup output directories and make sure pseudo-it will behave how you want it to. |
| `--maponly` | Only do one iteration and stop after read mapping. |
| `--noindels` | Set this to not incorporate indels into the final assembly. |
| `--diploid` | Set this use IUPAC ambiguity codes in the final FASTA file. |
| `--keepall` | By default, pseudo-it keeps only the final files for each step of each iteration (BAM, VCF, FASTA and their respective indices). Set this option to keep all intermediate files. While this is the best way to ensure your runs can be resumed with different settings this will result many large files being saved (total of ~1TB for a 30X genome and 4 iterations). |
| `--keeponlyfinal` | By default, pseudo-it keeps only the final files for each step of each iteration (BAM, VCF, FASTA and their respective indices). Set this option to keep these files ONLY for the final iteration. While this minimizes storage space required, you will be unable to resume this run. |
| `--overwrite` | Set this to overwrite existing files (as opposed to -resume which skips steps that have files already writen). |
| `--quiet` | Set this flag to prevent psuedo-it from reporting detailed information about each step. |
| `--version` | Simply print the version and exit. Can also be called as `-version`, `-v`, or `--v`. |
| `-h` | Print a help menu and exit. Can also be called as `--help`. |

## Pre-indexing your reference genome

Pseudo-it requires a reference genome in FASTA format. It also requires that this file be indexed in several ways. This indexing is kept separate from pseudo-it to ensure that the user has the most up-to-date indices for their file without risk of pseudo-it overwriting them.

Before you run pseudo it you MUST run the following commands on your `reference.fa` file:

1. Reference the file with samtools:

        samtools faidx <reference.fa>

2. Reference the file with BWA:

        bwa index <reference.fa>

3. Create a sequence dictionary with Picard:

        picard CreateSequenceDictionary R=<reference.fa> O=<reference.dict>

## Resource allocation

In addition to our own efforts to speed up the programs that pseudo-it runs, those programs themselves may also have the option to specify multiple processes/threads. This means that the allocation of processes/threads across programs can be confusing, and this section will try to explain how this is done.

1. First, multiple processes can be specified in the main pseudo-it interface with the `-p` option. This specifies the MAX number of processes to be used to run BWA and HaplotypeCaller. Minimally, each set of reads provided (`-se`, `-pe1` and `-pe2`, and/or `-pem`) would have 1 process for read mapping with BWA and each large scaffold/chromosome in the reference assembly would have 1 process for variant calling with HaplotypeCaller. 

2. BWA mem, the algorithm used to map reads by default, also allows for multiple threads with its internal `-t` option. Based on the number of read sets provided and the number of processes specified with `-p`, pseudo-it will automatically determine how many threads BWA mem should use.

    For example, if one set of paired end reads were supplied (`-pe1` and `-pe2`) and `-p` was set to 9, BWA mem would use all 9 threads for that single read set.

    If all three sets of reads are supplied (`-se`, `-pe1` and `-pe2`, and `-pem`) and `-p` was set to 9, BWA mem run all three sets simultaneously, each with 9 / 3 = 3 threads.

    Alternatively, you can specify the number of threads reserved for all BWA mem processes with the `-bwa-t` option. In the case of all three sets of reads being supplied with `-p 9` BUT with `-bwa-t 2`, each read set will be mapped simultaneously, but this time with only 2 threads a piece.

3. GATK's HaplotypeCaller also has a multi-thread option: `--native-pair-hmm-threads`. It is not clear how much this option effects runtime, but the default value with GATK is 4 which we leave as default for pseudo-it. 

    This means that if you run pseudo-it with a reference genome with 10 large scaffolds/chromosomes, you would ideally want to specify `-p 40` so that each chromosome gets its own GATK process that has access to 4 threads for the default `--native-pair-hmm-threads` option.

    Alternatively, with a reference genome with 10 large scaffolds/chromosomes, you could set pseudo-it's `gatk-t 1` to override the default `--native-pair-hmm-threads` option within GATK and you would only need `-p 10`.

    In the case of a fragmented assembly (i.e. 10 large scaffolds/chromosomes + 1000s of small scaffolds), you may want to specify one or two extra processes so the shorter scaffolds are not held up while the longer ones are being run.

## Example commands

### An optimized command for all 3 types of reads and a fragmented assembly with 10 large scaffolds + 1000s of small scaffolds

    python pseudo_it.py -ref [reference genome FASTA file] -se [single-end reads FASTQ file] -pe1 [paired-end reads FASTQ file 1] -pe2 [paired-end reads FASTQ file 2] -pem [merged paired-end reads FASTQ file] -i 4 -p 48 -o [desired output directory]

This command allocates 48 total processes for pseudo-it to use. It would run all three read sets through BWA mem with 48 / 3 = 16 threads each and would subsequently spawn 12 GATK HaplotypeCaller processes, each with the default `--native-pair-hmm-threads 4` option (12 * 4 = 48). The 10 long scaffolds would each have a separate processes while the shorter scaffolds would run have 2 extra processes allocated for them.

### An less resource-intensive command for all 3 types of reads and a fragmented assembly with 10 large scaffolds + 1000s of small scaffolds

    python pseudo_it.py -ref [reference genome FASTA file] -se [single-end reads FASTQ file] -pe1 [paired-end reads FASTQ file 1] -pe2 [paired-end reads FASTQ file 2] -pem [merged paired-end reads FASTQ file] -i 4 -p 12 -gatk-t 1 -o [desired output directory]

This command adds the `gatk-t 1` option to reduce the number of threads for HaplotypeCaller, thus reducing the total number of processes needed to `-p 12`. BWA mem would now run on each read set with only 4 threads.

# FAQs

**How many iterations should I perform?**

This depends on the sequence divergence of your sample relative to your reference and the type of data you have. For exome data, I found that three iterations performs well for samples with ~7.5 million years of divergence. If you expect more or have quickly evolving loci (noncoding, etc.), you might need more.