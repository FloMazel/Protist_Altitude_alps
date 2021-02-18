###############################################################################################################
###################################### ASV inference + F/R reads merg     #####################################
###############################################################################################################

rm(list=ls())

library(dada2)
library(tidyverse)
source("Scripts/0_Functions.R")

# Choose error rate for demultiplicaiton steps (1 mismatch authorized on the barcode)
n_core=4

# Define path to fastq files and outputs 
path_demultiplexed_runs="/Users/fmazel/Data/SOMETALP/processed_fastq/Subseted_Renamed_filtered_truncated_pooled_reads/"
#path_demultiplexed_runs="/Users/fmazel/Data/SOMETALP/processed_fastq/Subseted_Renamed_filtered_truncated_pooled_reads_1error/"

runs=list.files(path_demultiplexed_runs)

path_output="/Users/fmazel/Data/SOMETALP/processed_fastq/ASV_filtering_outputs/"
dir.create(path_output)

path_mergers_folder="/Users/fmazel/Data/SOMETALP/processed_fastq/ASV_filtering_outputs/Dada2_mergers/"
dir.create(path_mergers_folder)

for (irun in runs){
  # define path to run and list fastq files 
  pathFilt <- file.path(path_demultiplexed_runs,irun)
  
  filtFs <- grepl_wrap("_F_filt.fastq",list.files(pathFilt,full.names = T))
  filtRs <- grepl_wrap("_R_filt.fastq",list.files(pathFilt,full.names = T))
  
  #sampole names 
  sample.names <- sapply(strsplit(basename(filtFs), ".",1), `[`, 1)
  
  ####dereplication####
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  
  # Name the derep-class objects by the sample names
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  
  ####learn error rates####
  errF <- learnErrors(filtFs, multithread=n_core)
  errR <- learnErrors(filtRs, multithread=n_core)
  plotErrors(errF, nominalQ=TRUE)
  
  ####sample inference####
  dadaFs <- dada(derepFs, err=errF, multithread=n_core)
  dadaRs <- dada(derepRs, err=errR, multithread=n_core)
  
  ####merge paired reads####
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,minOverlap = 8) #a "subscript out of bounds" error here may indicate that you aren't merging any reads in one or more samples. you can remove samples with low counts from the workflow before the filterAndTrim step (a few steps back)
  
  saveRDS(mergers,paste(path_mergers_folder,irun,"_dada2_mergers.RDS",sep=""))
}


