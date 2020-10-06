#!/bin/bash
#SBATCH
#SBATCH --array=1-2%8
#SBATCH --partition=sixhour
#SBATCH --output=./slurmOutputs/seqTrim.%a.out
#SBATCH --error=./slurmOutputs/seqTrim.%a.out

# Script to trim primer sequences
# Assumes we are in the working directory
# Assumes working directory has the necessary subdirectories for storing output:
# ./trimmed/ # to hold trimmed sequence files
# ./slurmOutputs/ # to hold output

# Get input files (paired forward/reverse reads)
FWD=`ls ./raw/ | grep "_R1_" | head -n $SLURM_ARRAY_TASK_ID | tail -1` # get one forward read file
REV=`echo $FWD | sed 's/_R1_/_R2_/'` # get paired reverse read

# Extract sample name from the full input filename:
SAMPLE=$(echo ${FWD} | sed "s/_R1_\001\.fastq.gz*//")

# Fire up cutadapt
module load cutadapt/2.3

# Trim primers
### Forward primer = LROR = ACCCGCTGAACTTAAGC
### Reverse primer = FLR2 = TCGTTTAAAGCCATTACGTC
### Input comes from ./raw/ subdirectory
### Output goes into ./trimmed/ subdirectory
### Reads with no primer matches get discarded
cutadapt -g ACCCGCTGAACTTAAGC -G TCGTTTAAAGCCATTACGTC --discard-untrimmed --output=./trimmed/${SAMPLE}_R1_trimmed.fastq.gz --paired-output=./trimmed/${SAMPLE}_R2_trimmed.fastq.gz ./raw/$FWD ./raw/$REV
