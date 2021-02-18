#Analysis 


########################################################################
########## Within Groups dynamics for some parasite groups ############
########################################################################

rm(list=ls())

################################
####     Prep data          ####
################################

# Packages and functions
library(phyloseq)
library(tidyverse)
library(purrr)
library(cowplot)
library(RColorBrewer)
library(vegan)

source("Scripts/Final/Data_analysis/0_Functions.R")

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

cols=c('Oomycota'="darkblue",'Gregarinomorphea'="lightblue")
blues <- brewer.pal(8, "Blues")
cols=c('Oomycota'=blues[4],'Gregarinomorphea'=blues[8])


# Load data 
pooled_OTU_table = as.matrix(read.table("Data/Final/Pooled_Protists_ASVs_counts.txt",check.names = F))

pooled_metadata = read_csv("Data/Final/Pooled_metadata.csv") %>%
  mutate(plotID=as.character(plotID)) %>%
  mutate(altitude_bin=(cut_width(altitude, width=1000, boundary=0))) %>%
  mutate(pH_bin=cut_width(pH, width=1, boundary=7))


pooled_metadata$pH_bin

protists_taxonomy <- read_csv("Data/Final/Protists_taxonomy_functions.csv") %>%
  mutate(FUNCTIONS=ifelse(is.na(FUNCTIONS),'Unknown',FUNCTIONS))

alliance_display = expand_grid(color=brewer.pal(n = 5, name="Set1"),shape=c(15,17,18,25,8))[1:19,] 

# set up parameters

groups_taxo <- expand_grid(  rar=c(F,350),
                             taxo_rank= "PR2_Class",
                             taxo_group=c('Oomycota','Gregarinomorphea'),
                             minRead=400) %>% 
  mutate(file_output=ifelse(rar==F,"Redaction/V5/Figures/Fig_6.pdf","Redaction/V5/Figures/Supp_Figure_5.pdf"))

# Create phyloseq object 
MT <- as.data.frame(pooled_metadata); rownames(MT)=MT$plotID
PT <- as.matrix(protists_taxonomy); rownames(PT)=PT[,'seqs']
protist_PS=phyloseq(tax_table(PT),otu_table(pooled_OTU_table,taxa_are_rows=T),sample_data(MT))

# Info on alpha
protist_richness_non_rarefied <- map2_dfc(unique(groups_taxo$taxo_rank),unique(groups_taxo$taxo_group),estimate_richness_by_taxo_ranks,protist_PS) %>% 
  mutate(rarefied=NA) %>%
  mutate(plotID=sample_names(protist_PS)) %>%
  mutate(tot_depth = apply(as.matrix(otu_table(protist_PS)),2,sum))

# Add alpha on metadata
updated_metadata = pooled_metadata %>%
  left_join(protist_richness_non_rarefied)

summary(updated_metadata$readCounts_Oomycota)
summary(updated_metadata$Observed_Oomycota)

summary(updated_metadata$readCounts_Gregarinomorphea)
summary(updated_metadata$Observed_Gregarinomorphea)

MT <- as.data.frame(updated_metadata); rownames(MT)=MT$plotID
protist_PS=phyloseq(tax_table(PT),otu_table(pooled_OTU_table,taxa_are_rows=T),sample_data(MT))


################################
####    ANalysis           ####
################################

# Parameters of the analysis

dim1=function(x){dim(x)[1]}
source("Scripts/Final/Data_analysis/0_Functions.R")

# Run beta diversity computation, permanova, ordinations and plots
beta_df <- groups_taxo  %>%
  mutate(
    physeq = pmap(list(taxo_rank,taxo_group,minRead,rar),subset_physeqPR2_class,protist_PS), #subset associated metadata and filter 
    meta=map(physeq,subset_meta),
    beta = map2(rar,physeq,estimate_Beta), #compute beta diversity 
    permanova = map2(beta,meta,permanova_wrapper), # run permanovas, implies subseting plot/PCR replicate so all 
    nmds_2axes_model = map2(beta,meta,metaMDS_wrapper,k=2),
    nmds_2axes_stress = map(nmds_2axes_model,function(x){x$stress}),
    nmds_2axes_coord = map(nmds_2axes_model,extractMDS_coord),
    SampleSize= map(nmds_2axes_coord,dim1),# run NMDS
    all = map2(meta,nmds_2axes_coord,left_join,by='plotID'),
    altitude_plots = map2(all,taxo_group,plot_altitude_ordination),
    altitude_plots_legend = map2(all,taxo_group,plot_altitude_ordination_with_legend),
    alliance_plots = map2(all,taxo_group,plot_alliance_ordination),
    SoilWaterC_plots = map2(all,taxo_group,plot_soilWaterC_ordination))  

beta_df$physeq

stress = beta_df %>% 
  subset(file_output=="Redaction/V5/Figures/Fig_6.pdf") %>% 
  select(taxo_group,nmds_2axes_stress)

fig3output="Redaction/V5/Figures/Fig_6.pdf"

for (fig3output in unique(beta_df$file_output)){
  
  legend = get_legend(beta_df$altitude_plots_legend[[1]])# Extract legend
  
  plotlist <- beta_df %>% 
    subset(file_output==fig3output) %>% 
    pull(altitude_plots)
  #plotlist[[3]] <- legend
  
  # Plot NMDS altitud
  NMDS = plot_grid(plotlist =plotlist,
                   nrow =1,
                   labels = c("A. Oomycota","B. Gregarinomorphea"),
                   label_size = P3,
                   #labels = c("A","B", "C"),
                   hjust = c(-.1,-.1,-.1),
                    axis="b",align = "h")
  
  NMDS
  
  # plot Permanova reults 
  permanova_Results <- beta_df %>%
    subset(file_output==fig3output) %>% 
    select(taxo_group,permanova) %>%
    unnest(permanova) %>%
    mutate(Signif = ifelse(`Pr(>F)`>.05,"",ifelse(`Pr(>F)`>.01,"*",'**'))) %>% 
    mutate(predictorF = ifelse(predictor=="log(bulkSoilWaterContent)","Soil Water (%)",predictor)) %>% 
    mutate(predictorF = ifelse(predictor=="Alliance","Plants",predictorF)) %>% 
    mutate(predictorF = ifelse(predictor=="pH","Soil pH",predictorF)) %>% 
    mutate(predictorF = ifelse(predictor=="C.N","C/N ratio",predictorF)) %>% 
    mutate(predictorF = ifelse(predictor=="altitude","Elevation",predictorF)) 
  
  
  permanova_plot = permanova_Results %>% 
    subset(!is.na(F)) %>% 
    ggplot(aes(y=F,x=reorder(predictorF,desc(F)),fill=taxo_group))+
    geom_bar(stat='identity',position="dodge")+
    geom_text(aes(label=Signif), position=position_dodge(width=0.9), vjust=.5,size=6)+
    MyTheme+scale_fill_manual(values = cols)+#scale_color_manual(values = cols)+
    xlab('Environmental predictor')+ylab('Marginal Pseudo-F value')+
    xlab('Environmental predictor')+ylab('Marginal Pseudo-F value')+
    ggtitle(" ")+
    theme(legend.position="bottom")+guides(fill=guide_legend(title="Class",title.position = "top",title.hjust=.5))
  
  permanova_plot
  
  legendF = get_legend(permanova_plot)# Extract legend
  permanova_plot <- permanova_plot +
    theme(legend.position = "none") 
  
  permanova_plot
  
  R2_res = permanova_Results %>% 
    subset(predictor=="Residual") %>%
    mutate(R2=1-R2) 
  
  R2_plot = R2_res %>% 
    ggplot(aes(y=R2,x=taxo_group,fill=taxo_group))+
    geom_bar(stat='identity',position="dodge") +
    xlab('Class')+ylab('Model fit (r2)') + 
    MyTheme + ggtitle(" ")+scale_fill_manual(values = cols)
  
  R2_plot = R2_plot + theme(legend.position = "none")
  
  model_plot=plot_grid(permanova_plot,R2_plot,legend,nrow = 1,rel_widths = c(3,3,1),labels = c('D. Model predictors','E. Model fit'))
  
  model_plot=plot_grid(permanova_plot,R2_plot,
                       nrow = 1,
                       rel_widths =c(3,1),
                       labels = c('C','D'),
                       hjust = c(-.1,-.1),
                       axis="b",align = "h",
                       label_size = P3)
  
  model_plot
  
  legend_plot=plot_grid(legend,legendF,
                        rel_widths = c(1,1))
  legend_plot
  
  #Assemble figure 
  Fig_3 = plot_grid(NMDS,NULL,model_plot,NULL,legend_plot,nrow = 5,rel_heights = c(2,.1,2,.1,.5))
  Fig_3
  
  
  ggsave(filename = fig3output,device = 'pdf',width = 175,height = 175,units="mm")        
}



colOO = c(brewer.pal(12,"Set3"),brewer.pal(5,"Dark2"))
colApi = c(brewer.pal(12,"Set3"),brewer.pal(6,"Dark2"))


# Taxa plots


OOmTaxaBars = beta_df %>% 
  subset(taxo_group=="Oomycota"&rar==0) %>% 
  select(physeq) %>% 
  pull() %>% `[[`(1) %>% 
  transform_sample_counts(function(x) x / sum(x)) %>% 
  merge_samples("pH_bin") %>% 
  transform_sample_counts(function(x) x / sum(x)) %>% 
  psmelt() %>% 
  mutate(pH_bin=factor(Sample,levels=levels(pooled_metadata$pH_bin)),
         PR2_Genus=ifelse(is.na(PR2_Genus),"Unknowm",PR2_Genus)) %>% 
  ggplot(aes(y=Abundance,x=pH_bin,fill=PR2_Genus)) + 
  geom_bar(stat="identity", position = "stack") +
  scale_fill_manual(values=colOO) +
  MyTheme +ggtitle("B. Oomycota Genus composition along pH")+
  guides(fill=guide_legend(nrow=6))

OOmTaxaBars

Apicomplexa_taxaBars = beta_df %>% 
  subset(taxo_group=="Gregarinomorphea"&rar==0) %>% 
  select(physeq) %>% 
  pull() %>% `[[`(1) %>% 
  transform_sample_counts(function(x) x / sum(x)) %>% 
  merge_samples("altitude_bin") %>% 
  transform_sample_counts(function(x) x / sum(x)) %>% 
  psmelt() %>% 
  mutate(altitude_bin=factor(Sample,levels=levels(pooled_metadata$altitude_bin)),
         PR2_Genus=ifelse(is.na(PR2_Genus),"Unknowm",PR2_Genus)) %>% 
  ggplot(aes(y=Abundance,x=altitude_bin,fill=PR2_Genus)) + 
  geom_bar(stat="identity", position = "stack") +
  scale_fill_manual(values=colApi) +
  MyTheme +ggtitle("A. Gregarinomorphea Genus composition along altitude")+
  guides(fill=guide_legend(nrow=6))

Apicomplexa_taxaBars

taxa_p = plot_grid(Apicomplexa_taxaBars,
                   OOmTaxaBars,
                   nrow = 2)

taxa_p 
ggsave(plot = taxa_p ,filename = "Redaction/V5/Figures/Supp_Figure_6.pdf",device="pdf",height = 175,width = 175,units="mm")
