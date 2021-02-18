######################################################################################################################
###################################### OTU tables and fasta file produciton      #####################################
######################################################################################################################

rm(list=ls())

library(dada2)
library(tidyverse)
source("Scripts/0_Functions.R")

# Define paths
path_to_output="/Users/fmazel/Data/SOMETALP/processed_fastq/ASV_filtering_outputs/Final_outputs/"
dir.create(path_to_output)
source("Scripts/0_Functions.R")
mergers_folder_path="/Users/fmazel/Data/SOMETALP/processed_fastq/ASV_filtering_outputs/Dada2_mergers/"

# Import mergers 
Bacterial_mergers_path=list.files(mergers_folder_path,full.names = T)
Bacterial_mergers=unlist(lapply(Bacterial_mergers_path,readRDS),recursive = FALSE)

# Produce OTU table 
seqtab <- makeSequenceTable(Bacterial_mergers)

#look at sequence length distribution
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot

#chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=T, verbose=T) #
dim(seqtab.nochim)
dim(seqtab)

sum(seqtab.nochim)/sum(seqtab) # 10% of chimera, 

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {asv_headers[i] <- paste(">ASV", i, sep="_")}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste(path_to_output,"/Euk_ASVs.fa",sep=""))

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
#colnames(asv_tab)[1:50]=paste(colnames(asv_tab)[1:50],"replicate_1",sep="")
#colnames(asv_tab)[51:100]=paste(colnames(asv_tab)[51:100],"replicate_2",sep="")
#colnames(asv_tab)[101:150]=paste(colnames(asv_tab)[101:150],"replicate_3",sep="")
write.table(asv_tab, paste(path_to_output,"Euk_ASVs_counts.txt",sep=""), sep="\t", quote=F)

