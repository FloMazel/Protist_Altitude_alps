

##  Load data and code   ##
###########################
rm(list=ls())

#Load code
source("Scripts/Final/Data_analysis/0_Functions.R")
library(tidyverse)

#load data 
metadata = read_csv("Data/Final/Formated_Metadata_seq_outputs_allSites.csv") %>%
  mutate(altitude_bin=cut_width(altitude, width=500, boundary=0)) %>%
  mutate(plotID = as.character(plotID))

metadata_plot = read_csv("Data/Final/Formated_Metadata_byPlot_allSites.csv")

protists_taxonomy <- read_csv("Data/Final/Protists_taxonomy_functions.csv") %>%
  mutate(FUNCTIONS=ifelse(is.na(FUNCTIONS),'Unknown',FUNCTIONS))

  ASV_table=as.matrix(read.table("Data/Final/Protists_ASVs_counts.txt",check.names = F))

table(metadata$PCR_result_code)


########################################
#### CHOOSING THE SAMPLES TO KEEP ######
########################################

table(metadata$SampleWithinPlot)

subset_metadata = metadata %>%
  #mutate(sequencingID=gsub1(sequencingID,',','.')) %>%
  subset(PCR_result_code=="X") %>% # keep successful PCRs
  subset(!is.na(altitude)&SampleWithinPlot%in%c(NA,'.2','a',"T1")) # choice of sample when several sample per plot 

discarded_metadata = metadata %>%
  subset(PCR_result_code=="X") %>% # keep successful PCRs
  subset(!is.na(altitude)&!SampleWithinPlot%in%c(NA,'.2','a',"T1")) # choice of sample when several sample per plot 

nPCRrep = subset_metadata %>%
  group_by(plotID) %>%
  summarise(n_PCRrep=n())

metadata_plot = metadata_plot %>%
  left_join(nPCRrep)

colnames(ASV_table)[!colnames(ASV_table)%in%subset_metadata$UNIL_sequencingID]
subset_metadata$UNIL_sequencingID[!subset_metadata$UNIL_sequencingID%in%colnames(ASV_table)]

subset_OTU_table = ASV_table[,subset_metadata$UNIL_sequencingID]

########################################
####   Pool     PCR   replicates  ######
########################################

# Sites by sample matrix 
sitesampleM = subset_metadata %>%
  dplyr::select(plotID,UNIL_sequencingID) %>%
  mutate(X=1) %>%
  pivot_wider(id_cols = UNIL_sequencingID, values_from = X,names_from = plotID) %>% 
  tidyr::unnest()

sitesampleM 
sitesampleM1 = as.matrix(sitesampleM[,-1])
rownames(sitesampleM1) = sitesampleM$UNIL_sequencingID
sitesampleM1[is.na(sitesampleM1)]=0

# Pooling 
pooled_OTU_table = subset_OTU_table %*% sitesampleM1[colnames(subset_OTU_table),]
dim(pooled_OTU_table)

pooled_metadata = subset_metadata %>%
  group_by(plotID) %>%
  sample_n(1)

write.table(pooled_OTU_table,"Data/Final/Pooled_Protists_ASVs_counts.txt")
write_csv(pooled_metadata ,"Data/Final/Pooled_metadata.csv")




