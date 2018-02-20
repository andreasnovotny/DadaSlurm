#!/bin/Rscript -l

##########################################################################################################
############### DADA2 PIPELINE 2 : Learn Errors : Infer Sequence Variants ################################
##########################################################################################################
####                                                                                                  ####
#### Andreas Novotny, 2018-02                                                                         ####
#### https://github.com/andreasnovotny/DadaSlurm                                                                           ####
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
CURRENT_DIR <- args[1]
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


##########################################################################################################
#### 2. Learn forward and reverse error rates

print("R will now run the learnErrors function... ...")
errF <- learnErrors(filtFs, nread=2e6, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(filtRs, nread=2e6, multithread=TRUE, randomize=TRUE)

saveRDS(errF, file.path(CURRENT_DIR,'errF.rds'))
saveRDS(errR, file.path(CURRENT_DIR,'errR.rds'))

##########################################################################################################
#### 3. Sample inference and merger of paired-end reads

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

print("R will now loop derepFastq, dada and mergePairs functions... ...")
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR, maxMismatch=args[3], minOverlap=args[4])
    mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)


##########################################################################################################
#### 4. Construct sequence table
print("R will now makeSequenceTable... ...")
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file.path(CURRENT_DIR,'seqtab.rds'))


##########################################################################################################
##########################################################################################################
#### EXPLANATIONS: learnErrors

# Error rates are learned by alternating between sample inference and error rate estimation until convergence.
# Sample inferences is performed by the dada function. Error rate estimation is performed by
# errorEstimationFunction.
# The output of this function serves as input to the dada function call as the err parameter.
# Output: A named list with three entries:  $err_out:  A numeric matrix with the learned error rates.  $err_in:
# The initialization error rates (unimportant). $trans: A feature table of observed transitions for each
# type (eg. A->C) and quality score.

# OPTIONAL ARGUMENTS

# nreads
# Default: nreads = 1e+06. The minimum number of reads to use for error rate
# learning.  Samples are read into memory until at least this number of reads has
# been reached, or all provided samples have been read

# randomize
# Default FALSE. If FALSE, samples are read in the provided order
# until enough reads are obtained.  If TRUE, samples are picked at random from
# those provided

# MAX_CONSIST
# Default 10. The maximum number of times to step through the self-
# consistency loop. If convergence was not reached in MAX_CONSIST steps, the
# estimated error rates in the last step are returned.


##########################################################################################################
##########################################################################################################
#### EXPLANATIONS: derepFastq

# A custom interface to FastqStreamer for dereplicating amplicon sequences from fastq or com-
# pressed fastq files, while also controlling peak memory requirement to support large files.


##########################################################################################################
##########################################################################################################
#### EXPLANATIONS: dada

# The dada function takes as input dereplicated amplicon sequencing reads and returns the inferred
# composition of the sample (or samples).  Put another way, dada removes all sequencing errors to
# reveal the members of the sequenced community.
# If dada is run in selfConsist=TRUE mode, the algorithm will infer both the sample composition and
# the parameters of its error model from the data.s.

# Briefly, dada implements a statistical test for the notion that a specific sequence was seen too many
# times to have been caused by amplicon errors from currently inferred sample sequences.  Overly-
# abundant sequences are used as the seeds of new clusters of sequencing reads, and the final set of
# clusters is taken to represent the denoised composition of the sample. A more detailed explanation
# of the algorithm is found in two publications:
# •  Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP (2016).  DADA2:
#		High resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581-3.
# •  Rosen MJ, Callahan BJ, Fisher DS, Holmes SP (2012). Denoising PCR-amplified metagenome
#		data. BMC bioinformatics, 13(1), 283.

# dada depends on a parametric error model of substitutions. Thus the quality of its sample inference
# is affected by the accuracy of the estimated error rates.
# selfConsist mode allows these error rates to be inferred from the data.
# All comparisons between sequences performed by dada depend on pairwise alignments. This step
# is the most computationally intensive part of the algorithm, and two alignment heuristics have been
# implemented for speed:  A kmer-distance screen and banded Needleman-Wunsch alignmemt.  Sees etDadaOpt.

##########################################################################################################
##########################################################################################################
#### EXPLANATIONS: mergePairs

# This function attempts to merge each denoised pair of forward and reverse reads, rejecting any pairs
# which do not sufficiently overlap or which contain too many (>0 by default) mismatches in the
# overlap region.  Note:  This function assumes that the fastq files for the forward and reverse reads
# were in the same order.

# OPTIONAL ARGUMENTS

# minOverlap
# Default 20. The minimum length of the overlap required for merging
# the forward and reverse reads.

# maxMismatch
# Default 0. The maximum mismatches allowed in the overlap region.

# returnRejects
# Default FALSE. If TRUE, the pairs that that were rejected based on
# mismatches in the overlap region are retained in the return
# data.frame.

# justConcatenate
# Default FALSE. If TRUE, the forward and reverse-complemented
# reverse read are concatenated rather than merged, with a NNNNNNNNNN (10
# Ns) spacer inserted between them.

##########################################################################################################
##########################################################################################################
#### EXPLANATIONS: makeSequenceTable

# This function constructs a sequence table (analogous to an OTU table) from the provided list of samples.

# OPTIONAL ARGUMENTS

# orderBy
# (Optional). character(1).  Default "abundance".  Specifies how the sequences
# (columns) of the returned table should be ordered (decreasing).  Valid values:
# "abundance", "nsamples", NULL.
