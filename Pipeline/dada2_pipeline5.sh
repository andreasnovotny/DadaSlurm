#!/bin/bash -l

#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -o /proj/snic2017-7-248/nobackup/private/andreas/RD-1794/Outputs/slurm-%j.out
#SBATCH -J dada2_pipeline_5-8

##########################################################################################################
############### DADA2 PIPELINE 5-8: WORKFLOW FOR PAIR-END ILLUMINA AMPLICONS #############################
##########################################################################################################
####                                                                                                  ####
#### Andreas Novotny, 2018-02                                                                         ####
#### https://github.com/andreasnovotny/DadaSlurm                                                      ####
####                                                                                                  ####
##########################################################################################################
##########################################################################################################
##########################################################################################################

echo 'Running DadaSlurm pipeline 5-8'

##########################################################################################################
#### 6a. Make fasta file of seqtab sequences (CSVtoFASTA)
module load python3
python3 ${current_dir}/dada2_pipeline6.py $current_dir

##########################################################################################################
#### 6b. Align sequences with MUSCLE
module load bioinfo-tools
module load muscle
muscle -in ${current_dir}/seqs.fa -out ${current_dir}/seqs.afa

##########################################################################################################
#### 7. Build tree with PHANGORN
module load R_packages/3.4.3
Rscript ${current_dir}/dada2_pipeline7.R $current_dir

##########################################################################################################
#### 8. Ad tree to the phyloseq object
module load R_packages/3.4.3
Rscript ${current_dir}/dada2_pipeline8.R $current_dir

##########################################################################################################
##########################################################################################################
