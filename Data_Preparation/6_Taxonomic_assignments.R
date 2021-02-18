###############################################################################################################
######################################       Taxonomic assignments        #####################################
###############################################################################################################

rm(list=ls())

library(dada2)
library(tidyverse)
library(seqinr)
source("Scripts/0_Functions.R")
n_core=3 # number of core for parallel computing 

# Define path to ALL fasta files (fast)
ASV_fasta_files_path="/Users/fmazel/Data/SOMETALP/processed_fastq/ASV_filtering_outputs/Final_outputs/Euk_ASVs.fa"

# Load fasta 
seqs = unlist(read.fasta(ASV_fasta_files_path,as.string = TRUE, set.attributes = FALSE,forceDNAtolower=FALSE))

#### RDP within dada2 (slow) 

# assign taxonomy 
path_to_output="/Users/fmazel/Data/SOMETALP/processed_fastq/ASV_filtering_outputs/Final_outputs/"

taxo_PR2 <- assignTaxonomy(seqs = seqs, refFasta = '/Users/fmazel/Data/dada2_formated_databases/pr2_version_4.12.0_18S_dada2.fasta', taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"),outputBootstraps = T,multithread = T)
saveRDS(taxo_PR2,'Tmp_Taxonomy_RDP_PR2.RDS')
taxo_SILVA <- assignTaxonomy(seqs = seqs, refFasta = '/Users/fmazel/Data/dada2_formated_databases/silva_132.18s.99_rep_set.dada2.fa',outputBootstraps = T,multithread = T)

taxo_S_name = taxo_SILVA$tax; colnames(taxo_S_name) <- paste('SILVA_',colnames(taxo_S_name),sep="")
taxo_PR2_name = taxo_PR2$tax; colnames(taxo_PR2_name) <- paste('PR2_',colnames(taxo_PR2_name),sep="")

taxo_S_boot = taxo_SILVA$boot; colnames(taxo_S_boot) <- paste('Boot_SILVA_',colnames(taxo_S_boot),sep="")
taxo_PR2_boot = taxo_PR2$boot; colnames(taxo_PR2_boot) <- paste('Boot_PR2_',colnames(taxo_PR2_boot),sep="")

# create final data frame 
full_taxo = as_tibble(cbind(taxo_PR2_name,taxo_PR2_boot,taxo_S_name,taxo_S_boot))
full_taxo$seqs=seqs

write.csv(full_taxo,'/Users/fmazel/Data/SOMETALP/processed_fastq/ASV_filtering_outputs/Final_outputs/Taxonomy_RDP_PR2_SILVA.csv')

