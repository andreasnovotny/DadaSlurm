# DADA2-SLURM PIPELINE

BIOINFORMATIC ANALYSIS OF PAIR-END ILLUMINA AMPLICONS

https://github.com/andreasnovotny/DadaSlurm


## Description
This is a bioinformatic pipeline for automated and fast analysis of pair end Illumina data, starting with raw fastq sequences, ending with a full Phyloseq object m(McMurdie & Holmes 2013).

The first part of the pipeline, implementing the DADA2 R package (Callahan et al. 2016) is based on the tutorial: dada2 pipeline for big data by the same authors. The full tutorial and package description is available at: https://benjjneb.github.io/dada2/tutorial.html.

The second part of the pipeline uses the MUSCLE software to align the resulting sequences from the dada2 analysis (Edgar 2004), and Phangorn R package to construct a phylogenetic tree (Schliep 2011).


The pipeline is streamlined for use on the UPPMAX computational cluster Rackham.Visit https://www.uppmax.uu.se/ for more information.

#### Overview of the Pipeline components:

1. DADA2::filterAndTrim
2. DADA2::learnErrors, dada, mergePairs
3. DADA2::removeBimeraDenovo, assignTaxonomy
4. Phyloseq::phyloseq
5. (Bash linking script)
6. MUSCLE alignment
7. Phangorn::dml, NJ, fit
8. Phyloseq




## Installation
The pipeline needs no installation. Clone or download this github repository: https://github.andreasnovotny/DadaSlurm

Required software: R/3.4.3 with modules dada2/1.6 and phyloseq/1.22.3. The software and required packages are installed on the UPPMAX module system.



## Usage

#### 1. Input data
**Sequences:**
A file-path to demultiplexed, gzipped FASTQ sequences divided into a FWD and a REV directory.

**Database:**
A file-path to a learning database (see dada2 documentation).

**Metadata:**
A file-path to a .CSV file with row names corresponding to sequence sample names.

####  2. Modify the script
Use the master script: `dada2_pipeline.sh` do monitor the pipeline. Specify required pathways according to the instructions in the script. Several settings can be modified from the `dada2_pipeline.sh` directly and the pipeline can be run without modifying any of the other files.

Different components of the pipeline may be turned off using a `#` sign at the line of the component execution.

For detailed instructions regarding the dada2 pipeline, pleas see the package manual: https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf.

#### 3. Execute the pipeline
Copy the files of the pipeline into one directory
`scp PATH/TO/pipeline/dada2_pipeline* USERNAME@rackham.uppmax.uu.se:/proj/snic2017-x-xxx/nobackup/DIRECTORY` and execute the pipeline by writing `sbatch dada2_pipeline.sh`.

For more instructions on Working at UPPMAX, pleas read the file `Working_on_uppmax.md` in this repository.

## References

Benjamin J Callahan, Paul J McMurdie, Michael J Rosen, Andrew W Han, Amy Jo A Johnson & Susan P Holmes 2016. *DADA2: High-resolution sample inference from Illumina amplicon data*. Nature Methods volume 13, pages 581–583 (2016). doi:10.1038/nmeth.3869

Paul J. McMurdie, Susan Holmes 2013. *phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data* PlusOne, April 22, 2013. doi:10.1371/journal.pone.0061217

Edgar, R.C. (2004) *MUSCLE: a multiple sequence alignment method with reduced time and space complexity* BMC Bioinformatics, (5) 113


Klaus Peter Schliep 2011. *phangorn: phylogenetic analysis in R* ioinformatics, Volume 27, Issue 4, 15 February 2011, Pages 592–593. doi:10.1093/bioinformatics/btq706
