### *** README *** ###

##IMPORTANT NOTES:
#1. AMF LSU Pipeline
#2. primers : LROR (ACCCGCTGAACTTAAGC); FLR2 (TCGTTTAAAGCCATTACGTC)
#3. assumes demultiplexed files with names *R1_001.fastq.gz and *R2_001.fastq.gz
#4. Adapted from QIIME2 Tutorials and DADA2 pipeline modified from https://benjjneb.github.io/dada2/tutorial.html

###################################
####### Required programs #########
################################# be

# QIIME2/2019.10
# cutadapt/2.3
# fastqc/0.11.8
# R/3.6.1
# dada2

###################################
####### Directory structure #######
###################################

# Main working directory:
# could be named anything. All scripts use paths relative to this directory
mkdir myWD/ 
cd myWD/

# Subdirectory for raw sequence data:
mkdir ./raw/

# Subdirectory for primer-trimmed sequence data:
mkdir ./trimmed/

# Subdirectory for fastQC output files:
mkdir ./fastQCoutput

# Subdirectory for truncated + quality-filtered + concatenated sequence data:
mkdir ./filtered/

# Subdirectory for DADA2 output
mkdir ./dada2output

# Subdirectory for Qiime2 files
mkdir ./q2files

# Subdirectory for RAxML files
mkdir ./RAxMLfiles

# Subdirectory for all Output/Error files
mkdir ./slurmOutputs

###################################
#### Required files & locations ###
###################################

# *** All paths are relative to the working directory (.)

# Raw demultiplexed study sequences (paired-end):
./raw/*R1_001.fastq.gz
./raw/*R2_001.fastq.gz

# Reference database: 
./LRORFLR2.V10seqs_8.4.20.fasta

# Bash scripts to run on SLURM
./AMFseqTrimParallel.sh   # this one will trim primers from the raw sequences
./AMFtrimmedToOTUs.sh   # this one will quality-filter the primer-trimmed sequences, truncate reads to fixed lengths based on user input, denoise, concatenate FWD and REV reads into ASVs, and cluster into 97% OTUs.
./AMFbuildTree.sh   # this one places environmental (OTU) sequences onto the backbone tree using RAxML
./AMFlaunchTrees.sh   # this one optimizes AMFbuildTree.sh for parallel processing for a particular dataset (just decides how many parallel jobs need to be run for a given set of OTUs)

# R scripts:
./AMFdada2.R    # this gets called by AMFtrimmedToOTUs.sh
./AMFcladeExtract.R    # final script to run; outputs OTU table and rep-seqs file subsetted to include AMF only

###################################
######### Initial processing: #####
######### Remove primers and ######
######### visualize quality #######
###################################

nSamples=51
# IMPORTANT: update this to specify your number of samples

sbatch --array=1-$nSamples%8 AMFseqTrimParallel.sh		# Samples will be trimmed in parallel (up to 8 simultaneously)

### After AMFseqTrimParallel.sh is done running: 

# 1. Fire up fastQC:
module load fastqc

# 2. Generate fastQC summaries for R1 and R2:
zcat ./trimmed/*R1_trimmed.fastq.gz | fastqc -o ./fastQCoutput stdin:R1 
zcat ./trimmed/*R2_trimmed.fastq.gz | fastqc -o ./fastQCoutput stdin:R2

# 3. Download fastqc outputs (./fastQCoutput/*.html)
# 4. View on web browser: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# 5. Choose truncation lengths for R1 and for R2 -- 
#    where does sequence quality start to drop off?
#    Specify your cutoffs in the next sbatch command (below)

###################################
######### Quality filter, #########
######### denoise, and ############
######### make OTU table ##########
###################################

### *** NOTE: *** ###
## before proceeding, the following packages 
## must be installed in your R environment:
# Bioconductor/3.10
# dada2/1.15
# ShortRead
# Biostrings
# stringr

### Specify truncation lengths for R1 and R2:
# These (110 and 105) are default values; you should change them
# to fit your dataset (based on FastQC results)

sbatch --export=R1cutoff=110,R2cutoff=105 AMFtrimmedToOTUs.sh

# ^ this bash script will do everything from dada2 pipeline (make ASV table) to 
# clustering ASVs into 97% OTUs, to
# filtering the 97% OTU table in preparation for RAXML.

###################################
######### Place OTU sequences #####
###### on the backbone tree #######
###################################

# Place OTU representative sequences onto the backbone/reference tree using RAxML:
# This script will determine how many OTUs need to be placed, and launch a batch
# script with an appropriate number of tasks to be processed in parallel.

bash AMFlaunchTrees.sh

###################################
#### Separate AMF from non-AMF ####
#### OTUs and subset OTU table ####
###################################

### *** NOTE: *** ###
## before proceeding, the following packages 
## must be installed in your R environment:
# ape
# TreeTools

# After all tasks are finished running for AMFbuildTree.sh,
# run this R script to determine which OTUs fall within the
# AMF clade, and make subsets of the rep-seqs and OTU table
# for AMF-OTUs and non-AMF-OTUs

Rscript AMFcladeExtract.R
  