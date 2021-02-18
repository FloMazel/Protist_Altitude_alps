rm(list=ls())

#Load code
library(tidyverse)
library(ade4)
library("PerformanceAnalytics")

#Load data 
pooled_metadata = read_csv("Data/Final/Pooled_metadata.csv") %>%
  mutate(plotID=as.character(plotID))

# define variables 
edaphic <- c("soilTemp","bulkSoilWaterContent", "pH","C.N" ,"TOC","EC_1_5", "TotalP")
climatic <- c("bio1","bio10","bio11","bio12","bio13" ,"bio14", "bio15","bio16","bio17","bio18","bio19","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9")            

length(climatic)
length(edaphic)

# perform PCA 
factor_meta <- pooled_metadata %>% 
  select(altitude,all_of(edaphic),all_of(climatic)) %>% 
  drop_na()

PCA <- dudi.pca(factor_meta,nf = 3, scannf = FALSE)

pdf("Redaction/V4/Figures/Supp_Figure_1.pdf")
s.corcircle(PCA$co, lab = names(PCA$tab), full = FALSE, box = TRUE)
dev.off()

# Final choice of variables
edaphic <- c("bulkSoilWaterContent", "pH","C.N" ,"TOC")
climatic <- c("bio1")    
my_data <- factor_meta [, c('altitude',edaphic,climatic)]

pdf("Redaction/V4/Figures/Supp_Figure_2.pdf")
chart.Correlation(my_data, histogram=TRUE, pch=19,method = "spearman")
dev.off()

 





