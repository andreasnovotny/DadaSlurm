#!/bin/Rscript -l

##########################################################################################################
############################### DADA2 PIPELINE 1 : Filter and Trim #######################################
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

print("R is now running the filterAndTrim function... ...")
library(dada2); packageVersion("dada2")


##########################################################################################################
#### 1. File parsing
args <- commandArgs(TRUE)
OUTPUT_DIR <- args[1]
INPUT_DIR <- args[2]

pathF <- paste(INPUT_DIR,'/FWD', sep="")
pathR <- paste(INPUT_DIR,'/REV', sep="")

filtpathF <- file.path(pathF, "filtered")
filtpathR <- file.path(pathR, "filtered")

fastqFs <- sort(list.files(pathF, pattern=".fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern=".fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

##########################################################################################################
#### 2. Filter and Trim

#########################################################################################
#### MODIFY IMPORTANT PARAMETERS IN THE filterAndTrim FUNCTION!                     #####
#### https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf #####
#########################################################################################

filtered <- filterAndTrim(fwd=file.path(pathF, fastqFs),filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
							truncLen=c(270,250), #<------------------------------------------------------MODIFY!
							trimLeft=c(2,2), #<----------------------------------------------------------MODIFY!
							maxEE=c(1,1), #<-------------------------------------------------------------MODIFY!
							truncQ=0, #<-----------------------------------------------------------------MODIFY!
							maxN=0, rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)



saveRDS(filtered, file.path(OUTPUT_DIR,'/temporary/filtered.rds'))



##########################################################################################################
##########################################################################################################
