

##  Load data and code   ##
###########################
rm(list=ls())

#Load code
source("Scripts/Final/Data_analysis/0_Functions.R")
library(tidyverse)
library(phyloseq)
library(vegan)
library(cowplot)
library(RColorBrewer)

#define plot features 
# 1 pt = 0.35mm
# Recommended max size=7 , min size =5 in pt 
P1=5;P2=6;P3=7
MyTheme=theme_classic() + 
  theme(plot.title = element_text(size=P3),
        axis.text=element_text(size=P2),
        axis.text.x =element_text(size=P2,angle=45,hjust=1),
        axis.title=element_text(size=P2),
        
        legend.text=element_text(size=P2),
        legend.title=element_text(size=P2),
        
        panel.border = element_blank(),
        
        axis.ticks = element_line(size = 1*.35),
        axis.ticks.length = unit(.5, "mm"),
        
        plot.caption = element_text(hjust = 0), 
        plot.title.position = "plot")


cols=c('Phototrophic'="#FFD92F",'Parasite'="#1F78B4",'Consumer'="#A6761D",'Unknown'='grey')

# load data 
protists_taxonomy <- read_csv("Data/Final/Protists_taxonomy_functions.csv") %>%
  mutate(Functions=ifelse(is.na(Functions),'Unknown',Functions))

pooled_OTU_table = as.matrix(read.table("Data/Final/Pooled_Protists_ASVs_counts.txt",check.names = F))

pooled_metadata = read_csv("Data/Final/Pooled_metadata.csv") %>%
  mutate(plotID=as.character(plotID))

# Create phyloseq object 
# ------------------------

MT <- as.data.frame(pooled_metadata); rownames(MT)=MT$plotID
PT <- as.matrix(protists_taxonomy); rownames(PT)=PT[,'seqs']
FCT <- as.matrix(protists_taxonomy); rownames(FCT)=FCT[,'seqs']; FCT <- FCT[,c("Functions",'PR2_Class')]

protist_PS_FCT=phyloseq(tax_table(FCT),otu_table(pooled_OTU_table,taxa_are_rows=T),sample_data(MT))

FCT_PS <-  protist_PS_FCT %>%
  tax_glom(taxrank = 'Functions' ,NArm=F) %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>% 
  mutate(Functions=factor(Functions,levels = c("Unknown","Phototrophic","Parasite","Consumer")))


# Basic stats on local variability of FCT group abundance
FCT_PS %>% group_by(Sample,Functions) %>% 
  summarise(Ab=sum(Abundance)) %>% 
  subset(Functions == "Consumer") %>%
  pull(Ab) %>%
  summary()

FCT_PS %>% group_by(Sample,Functions) %>% 
  summarise(Ab=sum(Abundance)) %>% 
  subset(Functions == "Parasite") %>%
  pull(Ab) %>%
  summary()

FCT_PS %>% group_by(Sample,Functions) %>% 
  summarise(Ab=sum(Abundance)) %>% 
  subset(Functions == "Phototrophic") %>%
  pull(Ab) %>%
  summary()


# Plot the samples 
# ------------------------

minAb=0
fct_plot = FCT_PS %>%  filter(Abundance > minAb) %>% 
  mutate(Abu_Consumer=ifelse(Functions=='Consumer', Abundance,NA)) %>% 
  ggplot(aes(x = fct_reorder2(plotID, Functions, Abundance) , y = Abundance, fill =Functions)) + 
  geom_bar(stat = "identity") +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative read counts") + xlab("Samples (ordered by decreasing consumer abundance)") + scale_fill_manual(values=cols) +
  MyTheme + theme(axis.title.x = element_text(vjust=5),axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank()) 
fct_plot

ggsave("Redaction/V4/Figures/Fig_2.pdf",fct_plot,device = 'pdf', width = 175,height =60  ,units="mm")
