#!/bin/bash -l

#SBATCH -A snic2017-7-248
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 3:00:00
#SBATCH -J dada2_pipeline_1-4

##########################################################################################################
###############  DADA2 PIPELINE : WORKFLOW FOR PAIR-END ILLUMINA AMPLICONS ###############################
##########################################################################################################
####                                                                                                  ####
#### Andreas Novotny, 2018-02                                                                         ####
#### https://github.com/andreasnovotny/DadaSlurm                                                      ####
####                                                                                                  ####
#### Implemented from dada2 pipeline for big data                                                     ####
#### https://benjjneb.github.io/dada2/bigdata_paired.html                                             ####
####                                                                                                  ####
##########################################################################################################
##########################################################################################################
##########################################################################################################

# Copy to uppmax:     $ scp PATH/TO/dada2_pipeline*  USERNAME@rackham.uppmax.uu.se:/proj/snic2017-x-xxx/nobackup/DIRECTORY
# Start by using:     $ sbatch dada2_pipeline.sh
# Stop by using:      $ scancel JOBID
# Monitor by using:   $ jobinfo -u USERNAME
# Typically, the pipeline should be run on 16 cores for 6 h (Not including Muscle and Phangorn)

##########################################################################################################

# DEFINE YOUR ANALYSIS WORKING DIRECTORY
current_dir='/proj/snic2017-7-248/nobackup/private/andreas/P9607/testrun21'

#DEFINE YOUR RAW SEQUENCE DIRECTRY (containing required subdirectories /FWD and /REV)
sequence_dir='/proj/snic2017-7-248/nobackup/private/andreas/P9607/testseq'

#DEFINE YOUR DATABASE FILE.gz
database='/proj/snic2017-7-248/nobackup/private/andreas/P9607/pr2_version_4.72_dada2_Eukaryota.fasta.gz'
pr2_database=TRUE #Write TRUE if this is the pr2 database, otherwhise write FALSE

#DEFINE YOUR SAMPLE METADATA FILE.csv
metadata='/proj/snic2017-7-248/nobackup/private/andreas/P9607/metadata_18S.csv'


##########################################################################################################
##########################################################################################################
##########################################################################################################

module load R_packages/3.4.3

##########################################################################################################
#### 1. Filter sequences

#DEFINE filterAndTrim PARAMETERS:
Truncate_FWD=230 #Truncate  reads  after truncLen bases. Reads shorter than this are discarded
Truncate_REV=230
trimLeft_FWD=5 #The number of nucleotides to remove from the start of each read.
trimLeft_REV=5
maxEE_FWD=1 #After truncation, reads with higher than maxEE "expected errors" will be discarded. EE = sum(10^(-Q/10))
maxEE_REV=1
truncQ=5 #Truncate reads at the first instance of a quality score less than or equal to truncQ

Rscript ${current_dir}/dada2_pipeline1.R $current_dir $sequence_dir $Truncate_FWD $Truncate_REV $trimLeft_FWD $trimLeft_REV $maxEE_FWD $maxEE_REV $truncQ

##########################################################################################################
#### 2. Infer Sequence Variants

# DEFINE mergePairs PARAMETERS:
minOverlap=15 #The minimum length of the overlap required for merging the forward and reverse reads.
maxMismatch=1 #The maximum mismatches allowed in the overlap region.

Rscript ${current_dir}/dada2_pipeline2.R $current_dir $sequence_dir $maxMismatch $minOverlap

##########################################################################################################
#### 3. Remove Chimeras, Assign taxonomy

# DEFINE assignTaxonomy PARAMETERS:
minBoot=50 #The  minimum  bootstrap  confidence  for  assigning  a taxonomic level.

Rscript ${current_dir}/dada2_pipeline3.R $current_dir $database $pr2_database $minBoot

##########################################################################################################
#### 4. Combine data to a Phylosec object

Rscript ${current_dir}/dada2_pipeline4.R $current_dir $metadata

##########################################################################################################
#### 5. Optional: Construct Phylogeny With MUSCLE and PHANGORN

# DEFINE YOUR SLURM PROJECT ACCOUNT AGAIN:
A='snic2017-7-248'

sbatch -A $A --export=ALL,current_dir=$current_dir ${current_dir}/dada2_pipeline5.sh

##########################################################################################################

echo 'Finishing the dadaSlurm pipeline 1-4'

##########################################################################################################
##########################################################################################################
