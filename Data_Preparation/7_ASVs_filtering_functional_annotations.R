##################################################################
###  Clean ASVs data from non-Protists and assign Functions    ###
##################################################################

rm(list=ls()) 
library(tidyverse)
source('Scripts/Final/0_Functions.R')

# Load data
taxo <- read.csv('/Users/fmazel/Data/SOMETALP/processed_fastq/ASV_filtering_outputs/Final_outputs/Taxonomy_RDP_PR2_SILVA.csv',stringsAsFactors = F)
asv_table <- read.table("/Users/fmazel/Data/SOMETALP/processed_fastq/ASV_filtering_outputs/Final_outputs/Euk_ASVs_counts.txt")
row.names(asv_table)=as.character(taxo$seqs)

# Compute Total abundance 
ASV_abundance = apply(asv_table,1,sum)
ASV_prevalence = apply(asv_table>0,1,sum)
taxo$abundance = ASV_abundance
taxo$prevalence = ASV_prevalence

# Stats on unknown sequences
taxo %>% 
  subset(is.na(PR2_Division)) %>%
  pull(SILVA_Phylum) %>%
  table(useNA = 'always')

taxo %>%
  mutate(unknown_PR2_Division=is.na(PR2_Division)) %>%
  group_by(unknown_PR2_Division) %>%
  summarise(n=n(),
            totCounts=sum(abundance)) %>%
  mutate(relCounts=totCounts/sum(totCounts),
         relASVcounts=n/sum(n))

# Brief overview of non protist read abundance
taxo %>%
  mutate(non_protists=ifelse(PR2_Division%in%c('Fungi','Metazoa','Streptophyta'),PR2_Division,'Protist')) %>%
  group_by(non_protists) %>%
  summarise(n=n(),
            totCounts=sum(abundance)) %>%
  mutate(relCounts=totCounts/sum(totCounts),
         relASVcounts=n/sum(n))

#### ------------------------------------------------------------#### 
#### Remove non protist (fungi, metazoa, plants) and rare ASVs   #### 
#### ------------------------------------------------------------#### 

Protist_taxonomy = taxo %>% 
  subset(!PR2_Division%in%c('Fungi','Metazoa','Streptophyta') & # remove these 
        prevalence>1& #keep prevalents
        abundance>100) #keep abundants 

asv_table_protist = asv_table[as.character(taxo_protists$seqs),]
colnames(asv_table_protist) = sapply(colnames(asv_table_protist), extractseqnames)

write.table(asv_table_protist,'Data/new/Protists_ASVs_counts.txt')

#### ----------------------------------------------------#### 
###  Modify names for unknown taxonomy (replacing NA)    ####
#### ----------------------------------------------------#### 

for (i in 1: dim(Protist_taxonomy)[1]) {
  NAs=is.na(Protist_taxonomy[i,2:7])
  RankT = Protist_taxonomy[i,][max(c(2:7)[!NAs])]
  #Protist_taxonomy[i,2:7][NAs]
  if (sum(NAs>0)) {
    Protist_taxonomy[i,2:7][NAs]=paste('Unknown',RankT,sep="_")
  }
}

write_csv(Protist_taxonomy,file='Data/new/Protists_taxonomy.csv')

Protist_taxonomy = read_csv('Data/new/Protists_taxonomy.csv')

#### ----------------------------------------------------#### 
###            predict functions of ASVs                 ####
#### ----------------------------------------------------#### 

annotated_fct_groups  = read_csv("Results/Protists_annotation/Annotation_database/Final_annotation_keys.csv") %>% 
  select(Taxonomy,FUNCTIONS)

updated_Protist_taxonomy <- Protist_taxonomy %>%
  left_join(annotated_fct_groups, by=c('PR2_Class'='Taxonomy')) %>%  # Class level 
  left_join(annotated_fct_groups, by=c('PR2_Family'='Taxonomy')) %>%  # Family level 
  left_join(annotated_fct_groups, by=c('PR2_Genus'='Taxonomy')) %>%   # Genus level 
  mutate(Functions = ifelse(!is.na(FUNCTIONS),FUNCTIONS,ifelse(!is.na(FUNCTIONS.x),FUNCTIONS.x,ifelse(!is.na(FUNCTIONS.y),FUNCTIONS.y,NA)))) #%>% 

write_csv(updated_Protist_taxonomy,path='Data/Final/Protists_taxonomy_functions.csv')


