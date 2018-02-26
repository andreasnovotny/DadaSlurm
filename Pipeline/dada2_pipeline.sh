#!/bin/bash -l

#SBATCH -A snic201x-x-xxx
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 6:00:00
#SBATCH -J dada2_pipeline

##########################################################################################################
###############  DADA2 PIPELINE : ANALYSIS OF PAIR-END ILLUMINA AMPLICONS ###############################
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
current_dir='/proj/snic201X-x-xxx/nobackup/private/directory'

#DEFINE YOUR RAW SEQUENCE DIRECTRY (containing required subdirectories /FWD and /REV)
sequence_dir='/proj/snic201X-x-xxx/nobackup/private/directory'

#DEFINE YOUR DATABASE FILE.gz
database='/proj/snic201X-x-xxx/nobackup/private/database.fa.gz'
pr2_database=FALSE #Write TRUE if this is the pr2 database, otherwhise write FALSE

#DEFINE YOUR SAMPLE METADATA FILE.csv
metadata='/proj/snic201X-x-xxx/nobackup/private/metadata.csv'

# DEFINE YOUR SLURM PROJECT ACCOUNT AGAIN:
A='snic201X-x-xxx'


##########################################################################################################
##########################################################################################################
##########################################################################################################

module load R_packages/3.4.3

##########################################################################################################
#### 1. Filter sequences

Rscript ${current_dir}/dada2_pipeline1.R $current_dir $sequence_dir #$Truncate_FWD $Truncate_REV $trimLeft_FWD $trimLeft_REV $maxEE_FWD $maxEE_REV $truncQ

##########################################################################################################
#### 2. Infer Sequence Variants

Rscript ${current_dir}/dada2_pipeline2.R $current_dir $sequence_dir #$maxMismatch $minOverlap

##########################################################################################################
#### 3. Remove Chimeras, Assign taxonomy

Rscript ${current_dir}/dada2_pipeline3.R $current_dir $database $pr2_database #$minBoot

##########################################################################################################
#### 4. Create final output files.

Rscript ${current_dir}/dada2_pipeline4.R $current_dir $metadata

##########################################################################################################
#### 5. Optional: Construct Phylogeny With MUSCLE and PHANGORN

#sbatch -A $A --export=ALL,current_dir=$current_dir ${current_dir}/dada2_pipeline5.sh

##########################################################################################################

echo 'Finishing the dadaSlurm pipeline 1-4'

##########################################################################################################
##########################################################################################################
