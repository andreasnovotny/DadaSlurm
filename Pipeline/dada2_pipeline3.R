#!/bin/Rscript -l

##########################################################################################################
############### DADA2 PIPELINE 3 : Chimera Removal : Taxonomic Assignment ################################
##########################################################################################################
####                                                                                                  ####
#### Andreas Novotny, 2018-02                                                                         ####
#### https://github.com/andreasnovotny/DadaSlurm                                                                        ####
####                                                                                                  ####
#### Implemented from dada2 pipeline for big data                                                     ####
#### https://benjjneb.github.io/dada2/bigdata_paired.html                                             ####
####                                                                                                  ####
##########################################################################################################
##########################################################################################################
##########################################################################################################

library(dada2); packageVersion("dada2")

args <- commandArgs(TRUE)
CURRENT_DIR <- args[1]
DATABASE <- args[2]


##########################################################################################################
#### 1. File parsing

st.all <- readRDS(file.path(CURRENT_DIR,'seqtab.rds'))

##########################################################################################################
#### 2. Remove chimeras

print("R will now removeBimeraDenovo... ...")
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
saveRDS(seqtab, file.path(CURRENT_DIR,'seqtab_final.rds'))

##########################################################################################################
#### 3. Assign taxonomy

print("R will noe assignTaxonomy... ...")

if (args[3]==TRUE) {
  tax <- assignTaxonomy(seqtab, DATABASE, minBoot=args[4],
    outputBootstraps = TRUE,
    taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"),
    multithread=TRUE)
} else {
  tax <- assignTaxonomy(seqtab, DATABASE, minBoot=args[4],
    outputBootstraps = TRUE, multithread=TRUE)
}
saveRDS(tax, file.path(CURRENT_DIR,'tax_final.rds'))

##########################################################################################################
#### 4. Create big RSV table
seqtab_t <- t(seqtab)
RSV <- merge(tax, seqtab_t, by=0)
write.table(RSV, file.path(CURRENT_DIR,'RSV.txt'), sep = '\t')

##########################################################################################################
#### 5. make text file for MUSCLE alignment

seqs <- getSequences(seqtab)
write.csv(seqs, file.path(CURRENT_DIR,'seqs.csv'), sep = '\t')


##########################################################################################################
##########################################################################################################
#### EXPLANATIONS: removeBimeraDenovo

# This function is a convenience interface for chimera removal.  Two methods to identify chimeras
# are supported:  Identification from pooled sequences (see isBimeraDenovo for details) and identi-
# fication by consensus across samples (see isBimeraDenovoTable for details).  Sequence variants
# identified as bimeric are removed, and a bimera-free collection of unique sequences is returned.

# OPTIONAL ARGUMENTS

# method
# (Optional).  Default is "consensus".  Only has an effect if a sequence table is
# provided.
# If "pooled": The samples in the sequence table are all pooled together for bimera
# identification (isBimeraDenovo).
# If "consensus": The samples in a sequence table are independently checked for
# bimeras, and a consensus decision on each sequence variant is made (isBimeraDenovoTable).
# If "per-sample": The samples in a sequence table are independently checked for
# bimeras,  and sequence variants are removed (zeroed-out) from samples inde-
# pendently (isBimeraDenovo).

##########################################################################################################
##########################################################################################################
#### EXPLANATIONS: assignTaxonomy

# assignTaxonomy implements the RDP Naive Bayesian Classifier algorithm described in Wang et al.
# (Applied and Environmental Microbiology 2007), with kmer size 8 and 100 bootstrap replicates.
# Properly formatted reference files for several popular taxonomic databases are available
# http://benjjneb.github.io/dada2/training.html

# Output= A character matrix of assigned taxonomies exceeding the minBoot level of bootstrapping confi-
# dence. Rows correspond to the provided sequences, columns to the taxonomic levels. NA indicates
# that the sequence was not consistently classified at that level at the minBoot threshhold.
# If outputBootstraps is TRUE, a named list containing the assigned taxonomies (named "taxa") and
# the bootstrap values (named "boot") will be returned

# OPTIONAL ARGUMENTS

# minBoot
# Default  50.   The  minimum  bootstrap  confidence  for  assigning  a
# taxonomic level.

# tryRC
# Default FALSE. If TRUE, the reverse-complement of each sequences
# will be used for classification if it is a better match to the reference sequences
# than the forward sequence

# outputBootstraps
# Default FALSE. If TRUE, bootstrap values will be retained in an
# integer matrix. A named list containing the assigned taxonomies (named "taxa")
# and the bootstrap values (named "boot") will be returned.  Minimum bootstrap
# confidence filtering still takes place, to see full taxonomy set minBoot=0

# taxLevels
# (Optional).   Default  is  c("Kingdom",  "Phylum",  "Class",  "Order",  "Family",
# "Genus", "Species"). The taxonomic levels being assigned. Truncates if deeper
# levels not present in training fasta.
