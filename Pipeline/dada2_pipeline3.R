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

#########################################################################################
#### MODIFY IMPORTANT PARAMETERS IN THE assignTaxonomy FUNCTION!                    #####
#### https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf #####
#########################################################################################

print("R will now assignTaxonomy... ...")

if (args[3]==TRUE) {
  tax <- assignTaxonomy(seqtab, DATABASE, minBoot=50,
    outputBootstraps = TRUE,
    taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"),
    multithread=TRUE)
} else {
  tax <- assignTaxonomy(seqtab, DATABASE, minBoot=50,
    outputBootstraps = TRUE, multithread=TRUE)
}
saveRDS(tax, file.path(CURRENT_DIR,'tax_final.rds'))

##########################################################################################################
##########################################################################################################
