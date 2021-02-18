
##################################################
#####    Demultiplexing using cutadapt     #######
##################################################

rm(list=ls())

library(tidyverse)
library(seqinr)
library(readxl)

source("Scripts/0_Functions.R")

cutadapt <- "/Users/fmazel/Softwares/miniconda3/envs/cutadaptenv/bin/cutadapt" # path to cutadapt on the computer 
system2(cutadapt, args = "--version")

path_to_demultiplexed <- '/Users/fmazel/Data/SOMETALP/processed_fastq/Raw_demultiplexed/'
#path_to_demultiplexed <- '/Users/fmazel/Data/SOMETALP/processed_fastq/Raw_demultiplexed_1error/'
dir.create(path_to_demultiplexed)

path_to_metadata <- '/Users/fmazel/Data/SOMETALP/meta_physeq.csv'

path_to_raw_fastqfolder="/Users/fmazel/Data/SOMETALP/raw/PreAlps_Eukaryotes_UNINE/fastq_files/"
runs=list.files(path_to_raw_fastqfolder)

# --------------------------------#
# Create barcode fasta file in R  #
# --------------------------------#

barcodes=read.fasta("/Users/fmazel/Data/SOMETALP/raw/PreAlps_Eukaryotes_UNINE/V4_8nt_tags.fasta")
barcodes <- map(barcodes,function(x){c('^',toupper(x))}) # add > for cutadapt
path_to_barcode <- "/Users/fmazel/Data/SOMETALP/raw/PreAlps_Eukaryotes_UNINE/V4_8nt_tags_cutadapt.fasta"

barcode_names <- sub('-','_',names(barcodes))
barcode_names <- sub('V4','',barcode_names)

write.fasta(sequences=barcodes,names=barcode_names,file.out=path_to_barcode)
path_to_barcode <- paste('file:',path_to_barcode,sep='') 

# ----------------------------------#
# Demultiplex reads using cut adapt #
# ----------------------------------#

error=0.05

for (run in runs){
  print(run)
  
  fastq_files <- list.files(paste(path_to_raw_fastqfolder,run,sep=''),full.names = T) #path to fastq
  path_to_demultiplexed_run=paste(path_to_demultiplexed,run,sep="")
  dir.create(path_to_demultiplexed_run)
  
  params_cutadapt <- c('-e', error, # Allow no mismtaches in barcodes 
                       '--no-indels', # # No indels 
                       '-m', 100, # discard all reads < 100 pb
                       '-g', path_to_barcode, # Forward barcode 
                       '-G', path_to_barcode,  # Reverse barcode 
                       '-o',paste(path_to_demultiplexed_run,'/{name1}-{name2}.1.fastq.gz',sep=''), # R1 output
                       '-p',paste(path_to_demultiplexed_run,'/{name1}-{name2}.2.fastq.gz',sep=''), # R2 output 
                       fastq_files[1], # Input R1
                       fastq_files[2] # Input R2 
  )
  
  system2(cutadapt, args = params_cutadapt)
  
}
