# DADA2-SLURM PIPELINE

WORKFLOW FOR PAIR-END ILLUMINA AMPLICONS


## Description
Pipeline based on the tutorial: dada2 pipeline for big data https://benjjneb.github.io/dada2/bigdata_paired.html, streamlined for use on the UPPMAX cluster: https://www.uppmax.uu.se/.

The pipeline is a Slurm script implementing the following functions from the DADA2 Bioconductor R package.

1. filterAndTrim
2. learnErrors
3. derepFastq
4. dada
5. mergePairs
6. removeBimeraDenovo
7. assignTaxonomy

## Dependencies
Required software: R/3.4.3 with modules dada2/1.6 and phyloseq/1.22.3.
The software and required packages are installed on the UPPMAX module system.

## User Instructions
For detailed instructions regarding the dada2 pipeline, see the package manual: https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf.

#### 1. Installation
The pipeline needs no installation. Clone or download this github repository: https://github.andreasnovotny/DadaSlurm

#### 2. Input data
**Sequences:**
A file-path to demultiplexed, gzipped FASTQ sequences divided into a FWD and a REV directory.

**Database:**
A file-path to a learning database (see dada2 documentation).

**Metadata:**
A file-path to a .CSV file with row names coresponding to sequence sample names.

#### 3. Usage
The pipeline is monitored from the `pipeline.sh` bash script by following the instructions.

Typically, the pipeline should be run on 16 cores for 6 h.

Copy the pipeline to uppmax:
`scp PATH/TO/dada2_pipeline* USERNAME@rackham.uppmax.uu.se:/proj/snic2017-x-xxx/nobackup/DIRECTORY`

Start by using:     `sbatch dada2_pipeline.sh`

Stop by using:      `scancel JOBID`

Monitor by using:   `jobinfo -u USERNAME`



******************************************************

# Working on Uppmax - Instructions in summary.
Andreas Novotny
2018-02-19

andreas.novotny@su.se

### 1. Login

Login using the secure shell command:

`ssh *USERNAME*@rackham.uppmax.uu.se`

### 2. Upload and download files in uppmax

*Always run this commands from your lapdop terminal*

secure copy
`scp FROM: TO:`

Download file:
`scp *USERNAME*@rackham.uppmax.uu.se:/path/to/file.ext /local/path`

Uppload file:
`scp /local/path/to/file.ext *USERNAME*@rackham.uppmax.uu.se:/path/`

Directories:
Add the -r option:
`scp -r FROM: TO:`

### 3. The file system at Uppmax
`cd /path/to/dir/`

##### Your home directory

`cd ~`
or
`cd /home/*USERNAME*`

Only visible for the user. Dectory is backed up! Do NOT run any analysis here!

##### Project directories

`cd /proj/*PROJECT-NAME*`

-Directory is backed up! Do NOT run any analysis here - or you will run out of memory! Directory is public! Do NOT store nonpublic data here.

`cd /proj/*PROJECT-NAME*/private`

Dectory is backed up!Do NOT run any analysis here - or you will run out of memory! Directory is private! **USE THIS DIRECTORY FOR BACK UP STORAGE:**
	of raw data files; scripts; final final final results ect.

`cd /proj/*PROJECT-NAME*/nobackup`

Directory is public! Do NOT store nonpublic data here.

`cd /proj/*PROJECT-NAME*/nobackup/private`

Directory is only available for project members
Directory is NOT backed up **USE THIS DIRECORY FOR ANALYSIS:** temporary storage of analysis results ect.


### 4. The module system att uppmax:
Make modules available with:
`module load *MODULE*`

Important examples:

To Access all installed R packages:
`module load R_packages/3.4.3`

To Access all installed bioinformatics software:
`module load bioinfo-tools`
AND!
`module load *MODULE*`



### 5. The SLURM system for launcing big jobs:

Run analyses from Bash file scripts (myscript.sh)
Use the following template for SURM commands:

---
`#!/bin/bash -l`

`#SBATCH -A snic2017-x-xxx` <- the paying project

`#SBATCH -p core` <- core or node partition?

`#SBATCH -n 16` <- Number of cores

`#SBATCH -t 12:00:00` <- time (hh:mm:ss)

`#SBATCH -J dada2_pipeline` <- name the job

`module load python3` <- load needed modules

`FUNCTIONS.ect...`

---

Launch the script by using:
`sbatch /path/to/myscript.sh`

Stop process by:
`scancel *JOBID*`

Monitor jobs by:
`jobinfo -u *USERNAME*`
