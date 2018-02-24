#!/bin/Rscript -l

##########################################################################################################
############### DADA2 PIPELINE 4: Output Modifications ###################################################
##########################################################################################################
####                                                                                                  ####
#### Andreas Novotny, 2018-02                                                                         ####
#### https://github.com/andreasnovotny/DadaSlurm                                                      ####
####                                                                                                  ####
##########################################################################################################
##########################################################################################################
##########################################################################################################

library(phyloseq); packageVersion("phyloseq")
library(dada2)

args <- commandArgs(TRUE)
CURRENT_DIR <- args[1]
METADATA <- args[2]

print('R will now create a Phyloseq object from the results... ...')

##########################################################################################################
# 1. Read in all output files from the pipeline

seqtab_final <- as.matrix(readRDS(file.path(CURRENT_DIR,'seqtab_final.rds')))
taxonomy <- readRDS(file.path(CURRENT_DIR,'tax_final.rds'))
taxonomy1 <- as.matrix(taxonomy$tax)

sample.names <- readRDS(file.path(CURRENT_DIR, '/sample.names.rds'))
filtered <- readRDS(file.path(CURRENT_DIR,'/filtered.rds'))
mergers <- readRDS(file.path(CURRENT_DIR, '/mergers.rds'))
seqtab <- readRDS(file.path(CURRENT_DIR,'seqtab.rds'))

metadata <- read.csv(METADATA)
metadata2 <- metadata[,-1]
rownames(metadata2) <- metadata[,1]
metadata <- as.data.frame(metadata2)

##########################################################################################################
# 2. Create and save the phyloseq object

ps <- phyloseq(otu_table(seqtab_final, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxonomy1))
saveRDS(ps, file.path(CURRENT_DIR,'phyloseq.rds'))

##########################################################################################################
# 3. Track reads through the pipeline -make table:

track <- cbind(filtered, mergers, rowSums(seqtab), rowSums(seqtab_final))
colnames(track) <- c('input', 'filtered', 'mergd', 'tabled','nonchim')
rownames(track) <- sample.names

saveRDS(track, file.path(CURRENT_DIR,'/track_reads.rds'))

##########################################################################################################
#### 4. Create big RSV.csv table

seqtab_t <- t(seqtab_final)
RSV <- merge(taxonomy, seqtab_t, by=0)
write.csv(RSV, file.path(CURRENT_DIR,'RSV.txt'))

##########################################################################################################
#### 5. make csv file for MUSCLE alignment

seqs <- getSequences(seqtab_final)
write.csv(seqs, file.path(CURRENT_DIR,'seqs.csv'))

##########################################################################################################
##########################################################################################################
