#!/bin/bash -l

#SBATCH -A snic2017-7-248
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 2:00:00
#SBATCH -o /proj/snic2017-7-248/nobackup/private/andreas/RD-1794/Outputs/slurm-%j.out
#SBATCH -J dada2_pipeline

##########################################################################################################
###############  DADA2 PIPELINE : ANALYSIS OF PAIR-END ILLUMINA AMPLICONS ################################
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

# Start by using:     $ sbatch dada2_pipeline.sh
# Stop by using:      $ scancel JOBID
# Monitor by using:   $ jobinfo -u USERNAME
# Typically, the pipeline should be run on 16 cores for 6 h (Not including Muscle and Phangorn)

##########################################################################################################

# ANALYSIS NAME
analysis="18S_test"

# PATH TO RAW SEQUENCE DIRECTRY (containing required subdirectories /FWD and /REV)
sequence_dir='/proj/snic2017-7-248/nobackup/private/andreas/RD-1794/Sequences/testseq'

#PATH TO DATABASE FILE.gz
database='/proj/snic2017-7-248/nobackup/private/andreas/RD-1794/Metadata/pr2_version_4.72_dada2_Eukaryota.fasta.gz'
pr2_database=TRUE #Write TRUE if this is the pr2 database, otherwhise write FALSE

#PATH TO SAMPLE METADATA FILE.csv
metadata='/proj/snic2017-7-248/nobackup/private/andreas/RD-1794/Metadata/metadata.csv'

# PATH TO OUTPUT DIRECTORY
output_dir='/proj/snic2017-7-248/nobackup/private/andreas/RD-1794/Outputs'

# PATH TO PIPELINE DIRECTORY
pipeline_dir='/proj/snic2017-7-248/nobackup/private/andreas/RD-1794/DadaSlurm'

# PATH TO SLURM PROJECT ACCOUNT AGAIN
A='snic2017-7-248'


##########################################################################################################
##########################################################################################################
##########################################################################################################

module load R_packages/3.4.3

rm -r ${output_dir}/$analysis; mkdir ${output_dir}/$analysis
rm -r ${output_dir}/${analysis}/temporary; mkdir ${output_dir}/${analysis}/temporary
rm -r ${output_dir}/${analysis}/final; mkdir ${output_dir}/${analysis}/final
##########################################################################################################
#### 1. Filter sequences

Rscript ${pipeline_dir}/Pipeline/dada2_pipeline1.R ${output_dir}/$analysis $sequence_dir

##########################################################################################################
#### 2. Infer Sequence Variants

Rscript ${pipeline_dir}/Pipeline/dada2_pipeline2.R ${output_dir}/$analysis $sequence_dir

##########################################################################################################
#### 3. Remove Chimeras, Assign taxonomy

Rscript ${pipeline_dir}/Pipeline/dada2_pipeline3.R ${output_dir}/$analysis $database $pr2_database

##########################################################################################################
#### 4. Create final output files.

Rscript ${pipeline_dir}/Pipeline/dada2_pipeline4.R ${output_dir}/$analysis $metadata

##########################################################################################################
#### 5. Optional: Construct Phylogeny With MUSCLE and PHANGORN

#sbatch -A $A --export=ALL,current_dir=$current_dir ${pipeline_dir}/dada2_pipeline5.sh

##########################################################################################################

echo 'Finishing the dadaSlurm pipeline 1-4'

##########################################################################################################
##########################################################################################################
