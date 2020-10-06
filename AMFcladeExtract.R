#!/usr/bin/Rscript
####### Load libraries #######
library(ape)
library(TreeTools)

####### Determine which OTUs fall within the AMF clade #######
print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*') # make list of the tree files output from RAxML
for (tree in treeFiles) {
  print(paste0('Processing tree from OTU subset: ',tree))
  mytree <- read.tree(paste0('RAxMLfiles/',tree)) # load one tree
  # reformat tree to Preorder
  mytree.rooted <- RootTree(mytree,"EF417167.1_Rutilus_rutilus")
  mytree.rooted.preorder <- Preorder(mytree.rooted)
  # Make data frame with tip names and numbers
  tipDF <- data.frame('TipLabel'=TipLabels(mytree.rooted.preorder))
  tipDF$TipNumber <- rownames(tipDF) # copy row names (integer IDs) to column
  # get integer values of tips that define AMF clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='AM183920_Geosiphon_pyriformis')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832207_Acaulospora_tuberculata')$TipNumber)
  # Determine the node for the AMF clade:
  AMF.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF clade and identifying AMF OTUs:')
  # Extract the AMF clade:
  AMF.clade <- extract.clade(mytree.rooted.preorder,AMF.node)
  # plotTree(AMF.clade)
  # str(AMF.clade) # check how many tips
  ### Extract names of OTUs from the AMF clade:
  AMFnames<-c(AMF.clade$tip.label)
  ### Subset only those that start with 'OTU'
  hits<-AMFnames[grep("OTU",AMFnames)]
  if (tree==treeFiles[1]) { AMF.OTUs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMF.OTUs <- append(AMF.OTUs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF OTU names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMF.node, AMF.clade, AMFnames, hits) # clean up except for the AMF.OTUs vector
}

####### Load qiime2 output and prune out non-AMF OTUs #######
print('Loading qiime2 output to be pruned - OTU table')
# Load OTU table and prune:
otus <- read.delim('./otu97table_clean.tsv',sep='\t',skip=1,header=TRUE)
AMFotus <- subset(otus,X.OTU.ID %in% AMF.OTUs) # filter out non-AMF OTUs using the list we generated above
otherotus <- subset(otus,!(X.OTU.ID %in% AMF.OTUs)) # filter out non-AMF OTUs using the list we generated above
# Write to file:
write.table(AMFotus,file='./otu97table_clean_AMFonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)
write.table(otherotus,file='./otu97table_clean_nonAMF.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - OTU rep. seqs')
# Load OTU rep seqs and prune:
seqs <- read.FASTA('./otu97repseqs_clean.fasta', type = 'DNA')
AMFseqs <- seqs[AMF.OTUs]
otherseqs <- seqs[setdiff(names(seqs),AMF.OTUs)]
# Write to file:
write.FASTA(AMFseqs,'./otu97repseqs_clean_AMFonly.fasta')
write.FASTA(otherseqs,'./otu97repseqs_clean_nonAMF.fasta')

print('AMF - LSU pipeline complete.')
print('AMF-only OTU table and sequences have been written to file.')
print('non-AMF OTU table and sequences have also been written to file.')