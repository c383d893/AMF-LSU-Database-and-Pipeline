#!/usr/bin/Rscript

####### Load required packages #######
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(stringr)

####### File parsing #######

# set working directory:
wdPath <- getwd()
list.files(wdPath) # Make sure there are 'trimmed', 'filtered', and 'dada2output' subdirectories

# specify subdirectories
trimpath <- paste0(wdPath,"/trimmed/") # directory containing fwd and rev primer-trimmed reads
list.files(trimpath) # Make sure the trimmed files are in there

filtpath <- paste0(wdPath,"/filtered/") # directory containing fwd and rev quality-filtered reads
list.files(filtpath) # Should be empty at beginning of script

# This creates a list of file paths for the primer-trimmed input files
# It assumes forward and reverse fastq filenames have format:
# SAMPLENAME_R1_001.fastq.gz and SAMPLENAME_R2_001.fastq.gz
fastqFs.trimmed <- sort(list.files(trimpath, pattern="_R1_trimmed.fastq.gz", full.names = TRUE)) # get forward reads (R1)
fastqRs.trimmed <- sort(list.files(trimpath, pattern="_R2_trimmed.fastq.gz", full.names = TRUE)) # get reverse reads (R2)
if(length(fastqFs.trimmed) != length(fastqRs.trimmed)) stop("Forward and reverse files do not match.")

# Rewrite file names for filtered reads:
# It replaces 'trimmed' with 'filtered' both in the file name, and in the full path (subdirectory name)
fastqFs.filtered <- str_replace_all(fastqFs.trimmed,'trimmed','filtered')
fastqRs.filtered <- str_replace_all(fastqRs.trimmed,'trimmed','filtered')

# Extract simplified sample names
sample.names <- str_remove_all(basename(fastqFs.filtered),'_R1_filtered.fastq.gz')
names(fastqFs.filtered) <- sample.names
names(fastqRs.filtered) <- sample.names
print('Sample names are:\n')
print(sample.names)

####### Quality filtering + read truncation #######

# Set truncation lengths for fwd and rev reads
# The bash script will replace these values with user input
# Based on fastQC assessment
R1trunclen <- R1trunclen.value
R2trunclen <- R2trunclen.value

filterOutput <- filterAndTrim(fwd=fastqFs.trimmed, filt=fastqFs.filtered,
              rev=fastqRs.trimmed, filt.rev=fastqRs.filtered,truncLen=c(R1trunclen,R2trunclen),
              compress=TRUE, verbose=TRUE, multithread=TRUE) 
saveRDS(filterOutput,file=paste0(wdPath,'/dada2output/filterOutput.RDS'))

####### Learn the Error Rates and denoise #######
errF <- learnErrors(fastqFs.filtered, multithread=TRUE)
saveRDS(errF,file=paste0(wdPath,"/dada2output/errF.RDS")) # save intermediate file

errR <- learnErrors(fastqRs.filtered, multithread=TRUE)
saveRDS(errR,file=paste0(wdPath,"/dada2output/errR.RDS")) # save intermediate file

# Dereplicate:
derepFs <- derepFastq(fastqFs.filtered, verbose=TRUE) 
saveRDS(derepFs,file=paste0(wdPath,"/dada2output/derepFs.RDS")) # save intermediate file

derepRs <- derepFastq(fastqRs.filtered, verbose=TRUE)
saveRDS(derepRs,file=paste0(wdPath,"/dada2output/derepRs.RDS")) # save intermediate file

# Rename derep-class objects:
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Apply the core sample inference algorithm to the dereplicated data.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
saveRDS(dadaFs,file=paste0(wdPath,"/dada2output/dadaFs.RDS")) # save intermediate file

dadaRs <- dada(derepRs, err=errR, multithread=TRUE) 
saveRDS(dadaRs,file=paste0(wdPath,"/dada2output/dadaRs.RDS")) # save intermediate file

####### Merge paired-end reads by concatenating #######
# Merge the paired end reads, just concatenate
# instead of looking for overlapping regions.
# Will just stick both ends together with NNNNNNNNNN between.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate = TRUE, verbose=TRUE)
saveRDS(mergers,file=paste0(wdPath,"/dada2output/mergers.RDS"))

rm(derepFs,derepRs)   # clean up to save memory

#######  Construct sequence table ####### 
seqtab <- makeSequenceTable(mergers)

#######  Remove Chimeras ####### 
seqtab.nochim <- removeBimeraDenovo(
	seqtab, method="consensus", multithread=TRUE, verbose=TRUE
)

####### Report the number of reads that made it through each step in the pipeline: #######
getN <- function(x) sum(getUniques(x))
track <- cbind(
	filterOutput,
	sapply(dadaFs, getN),
	sapply(dadaRs, getN),
	sapply(mergers, getN),
	rowSums(seqtab.nochim)
)

colnames(track) <- c(
	"input",
	"filtered",
	"denoisedF",
	"denoisedR",
	"merged",
	"nonchim")
#rownames(track) <- sample.names
head(track)

####### write outputs #######

# Give sequence headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file=paste0(wdPath,"/dada2output/ASVs.fasta"))

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
head(asv_tab)

#add placeholder column, then row.names=FALSE
Taxa <- rownames(asv_tab)
head(Taxa)
data <- cbind(Taxa,asv_tab)
head(data)
write.table(
	data,
	file=paste0(wdPath,"/dada2output/ASVtable.tsv"),
	sep="\t", quote=F, row.names = FALSE
)

# Tracking table
write.table(
	track,
	file=paste0(wdPath,"/dada2output/tracking_table.tsv"),
	sep="\t", quote=F, col.names=NA
)