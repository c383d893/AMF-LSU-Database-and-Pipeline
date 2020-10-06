#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --partition=sixhour
#SBATCH --output=./slurmOutputs/AMFtrimmedToOTUs.out
#SBATCH --error=./slurmOutputs/AMFtrimmedToOTUs.out
#SBATCH --job-name=trm2otu

module purge  
module load R/3.6
module load qiime2
module list

# Replace placeholder text in R script with user-provided truncation lengths for R1 and R2
cat AMFdada2.R | sed "s/R1trunclen.value/$R1cutoff/" | sed "s/R2trunclen.value/$R2cutoff/" > AMFdada2withCutoffs.R

echo;echo "Beginning DADA2 pipeline..."
 Rscript AMFdada2withCutoffs.R
echo;echo "DADA2 pipeline complete. Converting output files to Qiime format..."
rm AMFdada2withCutoffs.R # remove temporary script once it's done running

# Convert ASV table from .tsv format to .biom format
biom convert -i ./dada2output/ASVtable.tsv -o ./q2files/ASVtable.biom --to-hdf5 --table-type="OTU table"
echo; echo "ASV table converted from .tsv to .biom."

# Convert ASV table from .biom to .qza, ASV sequences from .fasta to .qza
qiime tools import --input-path ./q2files/ASVtable.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ./q2files/ASVtable.qza
qiime tools import --input-path ./dada2output/ASVs.fasta --type 'FeatureData[Sequence]' --output-path ./q2files/ASVseqs.qza
# Convert reference database from .fasta to .qza
qiime tools import --input-path ./LRORFLR2.V10seqs_8.4.20.fasta --type 'FeatureData[Sequence]' --output-path ./q2files/AMFreferenceSeqs.qza
echo; echo "ASV table, ASV sequences, and reference sequences have been converted to .qza"

# OTU clustering
qiime vsearch cluster-features-open-reference --i-table ./q2files/ASVtable.qza --i-sequences ./q2files/ASVseqs.qza --i-reference-sequences ./q2files/AMFreferenceSeqs.qza --p-perc-identity 0.97 --o-clustered-table ./q2files/otu97table.qza --o-clustered-sequences ./q2files/otu97repseqs.qza --o-new-reference-sequences ./q2files/otu97newRefSeqs.qza --verbose
echo; echo "ASVs have been clustered into OTUs at 97% sequence similarity"

# OTU table cleanup 
# Remove OTUs that occur fewer than 5 times
qiime feature-table filter-features --i-table ./q2files/otu97table.qza --p-min-frequency 5 --o-filtered-table ./q2files/otu97table_minusRareOTUs_allSamples.qza --verbose
# Remove samples with fewer than 1000 reads
qiime feature-table filter-samples --i-table ./q2files/otu97table_minusRareOTUs_allSamples.qza --p-min-frequency 1000 --o-filtered-table ./q2files/otu97table_clean.qza --verbose
# Update OTU representative sequences to exclude OTUs that were removed:
qiime feature-table filter-seqs --i-data ./q2files/otu97repseqs.qza --i-table ./q2files/otu97table_clean.qza --o-filtered-data ./q2files/otu97repseqs_clean.qza --verbose
echo; echo "Rare OTUs and low-coverage samples have been removed from the dataset."

# Export OTU representative sequences into .fasta format:
# Qiime2 forces output into a file called dna-sequences.fasta . We direct it into ./q2files/
qiime tools export --input-path ./q2files/otu97repseqs_clean.qza --output-path ./q2files/
cat ./q2files/dna-sequences.fasta | sed 's/>/>OTU/g'|sed 's/ASV//g' > ./otu97repseqs_clean.fasta # rename file, relabel ASVs as OTUs,and any perfect matched to database append OTU and move to working directory.

# Export OTU table into .tsv format:
# Qiime2 forces output into a .biom file called feature-table.biom . We direct it into ./q2files/
qiime tools export --input-path ./q2files/otu97table_clean.qza --output-path ./q2files/
# Convert from .biom to .tsv:
biom convert -i ./q2files/feature-table.biom -o ./q2files/feature-table.tsv --to-tsv
cat ./q2files/feature-table.tsv | sed 's/ASV//g' | awk '{print "OTU"$0}' | sed 's/OTU#OTU ID/X.OTU.ID/g' > ./otu97table_clean.tsv # rename file, relabel ASVs as OTUs, any perfect matched to database append OTU and move to working directory                                              

echo; echo "AMFtrimmedToOTUs.sh has completed running. Be sure to check the out
put file for error messages."