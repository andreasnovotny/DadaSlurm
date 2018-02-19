#!/bin/Rscript -l

##########################################################################################################
############################### DADA2 PIPELINE 1 : Filter and Trim #######################################
##########################################################################################################
####                                                                                                  ####
#### Andreas Novotny, 2018-02                                                                         ####
#### andreas.novotny@su.se                                                                            ####
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
CURRENT_DIR <- args[1]
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

filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(args[3],args[4]), trimLeft=c(args[5],args[6]), maxEE=c(args[7],args[8]),
              truncQ=args[9], maxN=0, rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)



##########################################################################################################
##########################################################################################################
#### EXPLANATIONS: filterAndTrim

# Filters and trims an input fastq file(s) (can be compressed) based on several user-definable criteria,
# and outputs fastq file(s) (compressed by default) containing those trimmed reads which passed the
# filters.   Corresponding forward and reverse fastq file(s) can be provided as input,  in which case
# filtering is performed on the forward and reverse reads independently, and both reads must pass for
# the read pair to be output
# filterAndTrim is a multithreaded convenience interface for the fastqFilter
# and fastqPairedFilter filtering functions. Note that error messages and tracking are not handled gracefully when using the
# multithreading functionality.  If errors arise, it is recommended to re-run without multithreading to troubleshoot the issue.

# OPTIONAL ARGUMENTS

# truncLen
# Truncate  reads  after truncLen bases. Reads shorter than this are discarded.

# trimLeft
# The number of nucleotides to remove from the start of each read. If both truncLen and
# trimLeft are provided, filtered reads will have length truncLen-trimLeft.

# maxEE
# (Optional).  Default Inf (no EE filtering).  After truncation, reads with higher than
# maxEE "expected errors" will be discarded.  Expected errors are calculated
# from the nominal definition of the quality score: EE = sum(10^(-Q/10))

# truncQ
# Default 2. Truncate reads at the first instance of a quality score less than or equal to
# truncQ.

# maxN
# Default=0. How many N to allow. (dada2 does not allow any N)

# rm.phix
# Default TRUE. If TRUE, discard reads that match against the phiX genome, as determined by
# isPhiX.

# Compress
# makes output gzipped.

# verbose
# =TRUE to show output error mssage

# multithread
# Default is FALSE. If TRUE, input files are filtered in parallel via mclapply.  If an integer is provided, it is passed to the
# mc.cores argument of mclapply.  Note that the parallelization here is by forking, and each process is
# loading another fastq file into memory.  This option is ignored in Windows, as Windows does not support forking, with
# mc.cores set to 1.  If memory is an issue, execute in a clean environment and reduce the chunk size n
# and/or the number of threads.
