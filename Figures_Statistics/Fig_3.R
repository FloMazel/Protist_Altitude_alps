###############################
#   Protists communities      #
###############################

##  Load data and code   ##
###########################
rm(list=ls())

#Load code
library(tidyverse)
library(phyloseq)
library(vegan)
library(cowplot)
library(RColorBrewer)
source("Scripts/news_ASVs/0_Functions.R")

#define plot features 

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

minAb=0 
alliance_display = expand_grid(color=brewer.pal(n = 5, name="Set1"),shape=c(15,17,18,25,8))[1:19,] 

# load data 
protists_taxonomy <- read_csv("Data/Final/Protists_taxonomy_functions.csv") %>%
  mutate(Functions=ifelse(is.na(Functions),'Unknown',Functions))

pooled_OTU_table = as.matrix(read.table("Data/Final/Pooled_Protists_ASVs_counts.txt",check.names = F))

pooled_metadata = read_csv("Data/Final/Pooled_metadata.csv") %>%
  mutate(plotID=as.character(plotID),
         Alliance = gsub("_"," ",Alliance)) 

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


# Basic stats 
# ---------
alt_bin = FCT_PS %>% 
  mutate(altitude_bin1=cut_width(altitude, width=1000, boundary=0)) %>%
  group_by(altitude_bin1,Functions) %>% 
  summarise(mean_ab=mean(Abundance))

FCT_PS

# Figure 3 
# ---------

altitude_plot <- FCT_PS %>%  filter(Abundance > minAb) %>% 
  ggplot(aes(y=Abundance,x=altitude,color=Functions,shape=Functions))+geom_point(alpha=.2)+
  geom_smooth(method = "loess", fill = NA, show.legend = FALSE) + MyTheme + 
  scale_color_manual(values=cols) + scale_shape_manual(values=c(18,15, 16, 17))+
  ggtitle("A")+ theme(legend.position="none") +
  xlab("Elevation (meters ASL)")+ylab("Relative read counts") +   
  
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.box = "horizontal",
        axis.text.x =element_text(size=P2,angle=0,hjust=0.5)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=3)))

legend <- cowplot::get_legend(altitude_plot )

altitude_plot  <- altitude_plot  + theme(legend.position="none")
altitude_plot  



waterC_plot <- FCT_PS %>%  filter(Abundance > minAb) %>% 
  ggplot(aes(y=Abundance,x=bulkSoilWaterContent,color=Functions,shape=Functions))+geom_point(alpha=.4)+
  scale_x_continuous(trans='log10') +
  geom_smooth(method = "loess", fill = NA) +MyTheme +
  scale_color_manual(values=cols)+ scale_shape_manual(values=c(18,15, 16, 17))+
  ggtitle("B")+ ylab("") +xlab("Bulk soil water content (%)")+ theme(legend.position="none", axis.text.x =element_text(size=P2,angle=0,hjust=0.5))
waterC_plot

pH_plot <- FCT_PS %>%  filter(Abundance > minAb) %>% 
  ggplot(aes(y=Abundance,x=pH,color=Functions,shape=Functions))+geom_point(alpha=.2)+
  geom_smooth(method = "loess", fill = NA) +MyTheme+
  scale_color_manual(values=cols)+ scale_shape_manual(values=c(18,15, 16, 17))+
  ggtitle("C")+ theme(legend.position="none", axis.text.x =element_text(size=P2,angle=0,hjust=0.5))+
  ylab("Relative read counts")+xlab("soil pH")
pH_plot

alliance_fct_plot <- FCT_PS %>%  filter(!is.na(Alliance)) %>% 
  group_by(Alliance,Functions) %>%
  summarise(Abundance=mean(Abundance),altitude=mean(altitude)) %>%
  ggplot(aes(y=Abundance,x=fct_reorder2(Alliance,Functions, Abundance),fill=Functions))+#facet_wrap(~FUNCTIONS,ncol=1)+
  geom_bar(stat = "identity") +scale_fill_manual(values=cols)+
  MyTheme+
  xlab('Plant communities')+ ylab("")+ ggtitle("D")+
  theme(legend.title = element_blank(),
         legend.box = "horizontal",
         axis.ticks.x = element_blank(),
         axis.line.x = element_blank(),
         axis.line.y = element_blank())+ 
  guides(fill = guide_legend(nrow = 1,reverse = TRUE, keywidth = 1, keyheight = 1))+ 
  theme(legend.position="none")
  
alliance_fct_plot 


Upfct_env_plots <- plot_grid(altitude_plot,NULL,waterC_plot,
                           ncol=3, rel_widths=c(1,.1,1), align = "h")
lowfct_env_plots <- plot_grid(pH_plot,NULL,alliance_fct_plot,
                           ncol=3, rel_widths=c(1,.1,1))

fig3 <- plot_grid(Upfct_env_plots,NULL,lowfct_env_plots,NULL, legend,
                           ncol=1,rel_heights=c(1,.02,1,.02,.2))

fig3

ggsave("Redaction/V5/Figures/Fig_3.pdf",fig3,device = 'pdf',, width = 175,height =175  ,units="mm")


## Fig S1

Trait_Ratios_pooled <- FCT_PS %>%
  pivot_wider(id_cols = plotID, names_from = Functions,values_from = Abundance) %>%
  left_join(pooled_metadata,by="plotID") %>%
  mutate(Ratio_Parasites_Consumer=Parasite/Consumer,
          logRatio=-log(Ratio_Parasites_Consumer,10),
           Ratio_Consumer_parasite = 1/Ratio_Parasites_Consumer)

Data = Trait_Ratios_pooled %>% 
  subset(!is.na(Alliance)&Ratio_Consumer_parasite<50)  #%>%
#  group_by(Alliance) %>%
#  summarise(Alliance_SoilWaterC = median(SoilWaterC)) %>% ungroup()

# Get the residual from altitude
Fit = loess(log(1/Ratio_Parasites_Consumer) ~ altitude, data = Data, span = 0.75)
Data$residual_loess_altitude = Fit$residuals

library(scales)

Altitude_plot_with_alliances = Data %>%
  ggplot(aes(y=(1/Ratio_Parasites_Consumer),x=altitude,size=.5)) + 
  geom_point(aes(color=Alliance,shape=Alliance)) + 
  scale_shape_manual(values = alliance_display$shape) +
  scale_color_manual(values = alliance_display$color) +
  MyTheme + ylab("# Consumer reads / # Parasite reads") + xlab("Elevation (ASL)")+
  scale_y_continuous(trans='log',breaks = c(1,2,5,10,25,50))+
  ggtitle("A. Elevation") + guides(shape = guide_legend(nrow = 11),size= FALSE)
  Altitude_plot_with_alliances

Alliance_plot_residual_loessAltitute = Data %>%
  ggplot(aes(y=residual_loess_altitude,x=reorder(Alliance,desc(residual_loess_altitude)),color=log(bulkSoilWaterContent))) + 
  geom_boxplot(outlier.size = 0) + geom_jitter(size=3) +
  MyTheme + ylab("# Consumer reads / # Parasite reads \n(Residual from a loess fit with altitude" ) + xlab('Plant alliances')+
  scale_color_gradient(low = "gold", high = "blue", space = "Lab" )+
  ggtitle("B. Plant communities and Soil water Content (%)") + xlab("Plant communities")+
  guides(col=guide_legend(title="Soil W. Content (%, log scale)"))

Altitude_plot_with_alliances
Alliance_plot_residual_loessAltitute
fct_env_plots2 <- plot_grid(Altitude_plot_with_alliances ,Alliance_plot_residual_loessAltitute,ncol=1)

ggsave("Redaction/V4/Figures/Supp_Figure_3.pdf",fct_env_plots2,device = 'pdf', width = 188,height = 188,units="mm")
