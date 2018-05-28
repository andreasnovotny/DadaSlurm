#!/bin/Rscript -l

##########################################################################################################
############### DADA2 PIPELINE 2 : Learn Errors : Infer Sequence Variants ################################
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

library(dada2); packageVersion("dada2")

##########################################################################################################
#### 1. File parsing

args <- commandArgs(TRUE)
OUTPUT_DIR <- args[1]
INPUT_DIR <- args[2]

filtpathF <- paste(INPUT_DIR,'/FWD/filtered', sep="")
filtpathR <- paste(INPUT_DIR,'/REV/filtered', sep="")

filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_L001_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_L001_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

saveRDS(sample.names, file.path(OUTPUT_DIR, 'temporary/sample.names.rds'))

##########################################################################################################
#### 2. Learn forward and reverse error rates

#########################################################################################
#### MODIFY IMPORTANT PARAMETERS IN THE learnErrors FUNCTION!                       #####
#### https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf #####
#########################################################################################

print("R will now run the learnErrors function... ...")
errF <- learnErrors(filtFs,
	nread=2e6, #<-----------------------------------------------------------MODIFY!
	multithread=TRUE, randomize=TRUE)
errR <- learnErrors(filtRs,
	nread=2e6, #<-----------------------------------------------------------MODIFY!
	multithread=TRUE, randomize=TRUE)

saveRDS(errF, file.path(OUTPUT_DIR,'temporary/errF.rds'))
saveRDS(errR, file.path(OUTPUT_DIR,'temporary/errR.rds'))

##########################################################################################################
#### 3. Sample inference and merger of paired-end reads

#########################################################################################
#### MODIFY IMPORTANT PARAMETERS IN THE mergePairs FUNCTION!                        #####
#### https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf #####
#########################################################################################

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

print("R will now loop derepFastq, dada and mergePairs functions... ...")
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR,
			maxMismatch=1, #<---------------------------------------------------MODIFY!
			minOverlap=15) #<---------------------------------------------------MODIFY!
    mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

getN <- function(x) sum(getUniques(x))
mergetab <- sapply(mergers, getN)
saveRDS(mergetab, file.path(OUTPUT_DIR, 'temporary/mergers.rds'))





##########################################################################################################
#### 4. Construct sequence table
print("R will now makeSequenceTable... ...")
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file.path(OUTPUT_DIR,'temporary/seqtab.rds'))








##########################################################################################################
##########################################################################################################
