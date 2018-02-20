#!/bin/Rscript -l

##########################################################################################################
############### DADA2 PIPELINE 7 : Build Phylogeny with Phangorn #########################################
##########################################################################################################
####                                                                                                  ####
#### Andreas Novotny, 2018-02                                                                         ####
#### https://github.com/andreasnovotny/DadaSlurm                                                      ####
####                                                                                                  ####
##########################################################################################################
##########################################################################################################
##########################################################################################################

print("Running R PHANGORN script")
library(phangorn)

##########################################################################################################
# 1. File parsing

args <- commandArgs(TRUE)
CURRENT_DIR <- args[1]
input_alignment <- file.path(CURRENT_DIR,'seqs.afa')

##########################################################################################################
# 2.Convert .afa aligned fasta fro phyDat (phangorn) object
phang.align <- read.phyDat(file=input_alignment, format="fasta", type="DNA")

##########################################################################################################
# 3.Create distance matrix
print('PHANGORN will now create a distande matrix of the aligne sequences... ...')
dm <- dist.ml(phang.align)
saveRDS(dm, file_path(CURRENT_DIR,'dm.rds'))

##########################################################################################################
# 4.Neighbour joining tree
print('PHANGORN will now create a Neighbour joining tree from the distance matrix... ...')
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

saveRDS(treeNJ, file.path(CURRENT_DIR,'treeNJ.rds'))
saveRDS(fit, file.path(CURRENT_DIR,'fit.rds'))

##########################################################################################################
# 5. Generalized time-reversible with Gamma rate variation maximum likelihood tree
print('PHANGORN will now generate the final GTR tree... ...')
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
rearrangement = "stochastic", control = pml.control(trace = 0))
saveRDS(fitGTR, file.path(CURRENT_DIR,'phangorn_tree.rds'))

##########################################################################################################
##########################################################################################################
#### EXPLANATIONS:

# The phangorn R package is used to construct a phylogenetic tree.
# We first construct a neighbor-joining tree, and then fit a GTR+G+I
# (Generalized time-reversible with Gamma rate variation)
# maximum likelihood tree using the neighbor-joining tree as a starting point.
