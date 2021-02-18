###############################################################################################
###############       Select/Filter/Trim/Rename/Pool demultiplexed outputs      ############### 
###############################################################################################

rm(list=ls())

library(tidyverse)
library(dada2)
library(ShortRead)
library(Biostrings)

source("Scripts/0_Functions.R")

n_core = 4 

cutadapt <- "/Users/fmazel/Softwares/miniconda3/envs/cutadaptenv/bin/cutadapt" # path to cutadapt on the computer 
system2(cutadapt, args = "--version")

# Load Metadata
metadata <- read_csv("Data/Formated_Metadata.csv") %>%
  mutate(Forward_barcode=paste('F',toupper(`P1 F`),sep='_'), #change names of barcode code
         Reverse_barcode=paste('R',toupper(`P2 R`),sep='_'),  #change names of barcode code
         Barcode_combn = paste(Forward_barcode,Reverse_barcode,sep="__"), #combine pairs into one single name 
         Barcode_combn_2 = paste(Reverse_barcode,Forward_barcode,sep="__"),  #combine other orientation pairs into one single name 
         Plate=as.character(PCR_replicate)) %>% #add plate number 
         pivot_longer(cols=Barcode_combn:Barcode_combn_2,names_to = 'Barcode_orientation',values_to = 'Barcode_combo') # provide both sequence of barcodes (i.e. illumina reads 1 can contain F or R barcodes, same for illumina read 2)


# Define input and outputs directories 
path_demultiplexed_runs="/Users/fmazel/Data/SOMETALP/processed_fastq/Raw_demultiplexed/"
#path_demultiplexed_runs="/Users/fmazel/Data/SOMETALP/processed_fastq/Raw_demultiplexed_1error/"
list_demultiplexed_runs=list.files(path_demultiplexed_runs)

# Define path to folder with renamed and subseted fastq files 
path_to_subseted_cleaned_runs="/Users/fmazel/Data/SOMETALP/processed_fastq/Subseted_Renamed_filtered_truncated_reads/"
#path_to_subseted_cleaned_runs="/Users/fmazel/Data/SOMETALP/processed_fastq/Subseted_Renamed_filtered_truncated_reads_1error/"
dir.create(path_to_subseted_cleaned_runs)

# combine fastq (pool mixed orientation reads)
path_to_subseted_cleaned_pooled_runs="/Users/fmazel/Data/SOMETALP/processed_fastq/Subseted_Renamed_filtered_truncated_pooled_reads/"
#path_to_subseted_cleaned_pooled_runs="/Users/fmazel/Data/SOMETALP/processed_fastq/Subseted_Renamed_filtered_truncated_pooled_reads_1error/"
dir.create(path_to_subseted_cleaned_pooled_runs)

# Get statistics on demultiplexed samples
fastq_file_characteristics_raw <- map_dfr(list_demultiplexed_runs,import_analyze_runs,path_demultiplexed_runs,splitName=F)

fastq_file_characteristics = fastq_file_characteristics_raw %>% # Assign sampleIDs, sample description to demultiplexing outputs   
  separate(run,sep="_",into=c(NA,NA,'Date',NA,NA,'Plate',NA),remove = F)  %>%
  separate(Plate,sep="e",into=c(NA,'Plate')  ) %>%
  separate(fastqfile_name,"-",into=c('Forward_barcode','tmp'),remove=F) %>%
  separate(tmp,"([.])",into=c('Reverse_barcode','Read'),extra="drop") %>%
  mutate(Barcode_combn = paste(Forward_barcode,Reverse_barcode,sep="__")) %>%
  full_join(metadata,by=c('Barcode_combn'='Barcode_combo','Plate'='Plate')) %>% #Add samples decription 
  mutate(Barcode_code=ifelse(is.na(plotID),"Non_used_barcode",'Used_barcode'))

fastq_file_characteristics %>% 
  group_by(Barcode_code) %>%
  summarise(rel_counts=sum(read_number)/sum(fastq_file_characteristics_raw$read_number),
            tot_counts=sum(read_number)) # 20 millions reads per run

# Bacteria reads 
for (irun in list_demultiplexed_runs){ # For each run
  
  # Define path to diferent input/output files 
  path_to_fastq=paste(path_demultiplexed_runs,irun,sep="") #input 
  list_fastq=list.files(path_to_fastq) # List all demultiplexed files
  
  path_to_filt_fastq=paste(path_to_subseted_cleaned_runs,irun,sep="") # initial output after filtering 
  dir.create(path_to_filt_fastq)
  
  path_to_combn_fastq=paste(path_to_subseted_cleaned_pooled_runs,irun,"/",sep="") # final output after combining 
  dir.create(path_to_combn_fastq)
  
  # Subset to sometalp sample 
  subset_files <- fastq_file_characteristics %>%
    subset(run==irun&PCR_replicate%in%c(1,2,3))
  
  # Define corresponding fastq file paths
  Fs = subset(subset_files,Read==1)$fastqfile_name
  Rs = subset(subset_files,Read==2)$fastqfile_name
  fnFs <- paste(path_to_fastq,Fs,sep="/")
  fnRs <- paste(path_to_fastq,Rs,sep="/")
  
  #plot quality profiles for forward and reverse reads ####
  plotQualityProfile(fnFs[1:20])
  plotQualityProfile(fnRs[1:20])
  
  #trim & filter####
  filtFs <- file.path(path_to_filt_fastq, Fs)
  filtRs <- file.path(path_to_filt_fastq, Rs)
  
  #out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,200),
  #                     maxN=0, maxEE=c(6,6), truncQ=2, rm.phix=TRUE,matchIDs=F,
  #                     compress=TRUE, multithread=n_core) #the more aggressive you are with truncation, the more reads will be retained, but you still need to merge, so consider how much read overlap remains after truncation
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(0,0),
                       maxN=0, maxEE=c(6,6), truncQ=2, rm.phix=TRUE,matchIDs=F,
                       compress=TRUE, multithread=n_core) #the more aggressive you are with truncation, the more reads will be retained, but you still need to merge, so consider how much read overlap remains after truncation
  
  percent_retained <- as.data.frame(out)
  percent_retained$percentage <- (percent_retained[,2]/percent_retained[,1])*100
  hist(percent_retained$percentage)
  
  # Combine fastq files 
for (i in unique(subset_files$sequencingID)){
    sets <- subset(subset_files,sequencingID==i)  
    P5_F_1 <- paste(path_to_filt_fastq,subset(sets,grepl("F_",Forward_barcode.x)==T&Read==1)$fastqfile_name,sep="/")
    P5_F_2 <- paste(path_to_filt_fastq,subset(sets,grepl("F_",Forward_barcode.x)==T&Read==2)$fastqfile_name,sep="/")
    P7_F_1 <- paste(path_to_filt_fastq,subset(sets,grepl("R_",Forward_barcode.x)==T&Read==1)$fastqfile_name,sep="/")
    P7_F_2 <- paste(path_to_filt_fastq,subset(sets,grepl("R_",Forward_barcode.x)==T&Read==2)$fastqfile_name,sep="/")
    
    newP5_F_1 <- append(readFastq(P5_F_1),readFastq(P7_F_2)) # adding the reverse reads of the file with 806R in forward reads in new forwards reads
    newP5_F_2 <- append(readFastq(P5_F_2),readFastq(P7_F_1)) # adding the forward reads of the file with 806R in forward reads in new reverse reads
    
    writeFastq(newP5_F_1,paste(path_to_combn_fastq,i,"_F_filt.fastq",sep=""),mode='w')
    writeFastq(newP5_F_2,paste(path_to_combn_fastq,i,"_R_filt.fastq",sep=""),mode='w')
  }
  
  
}
