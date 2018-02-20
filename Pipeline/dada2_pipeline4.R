#!/bin/Rscript -l

##########################################################################################################
############### DADA2 PIPELINE 4: Create Phyloseq Object #################################################
##########################################################################################################
####                                                                                                  ####
#### Andreas Novotny, 2018-02                                                                         ####
#### https://github.com/andreasnovotny/DadaSlurm                                                      ####
####                                                                                                  ####
##########################################################################################################
##########################################################################################################
##########################################################################################################

library(phyloseq); packageVersion("phyloseq")

args <- commandArgs(TRUE)
CURRENT_DIR <- args[1]
METADATA <- args[2]

print('R will now create a Phyloseq object from the results... ...')

##########################################################################################################
# 1. Read in all output files from the pipeline

seqtab <- as.matrix(readRDS(file.path(CURRENT_DIR,'seqtab_final.rds')))
taxonomy <- readRDS(file.path(CURRENT_DIR,'tax_final.rds'))
taxonomy <- as.matrix(taxonomy$tax)

metadata <- read.csv2(METADATA)
metadata2 <- metadata[,-1]
rownames(metadata2) <- metadata[,1]
metadata <- as.data.frame(metadata2)

##########################################################################################################
# 2. Create and save the phyloseq object

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxonomy))

saveRDS(ps, file.path(CURRENT_DIR,'phyloseq.rds'))

##########################################################################################################
##########################################################################################################
