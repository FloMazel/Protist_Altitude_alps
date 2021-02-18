
##  Load data and code   ##
###########################
rm(list=ls())

#Load code
library(tidyverse)
library(phyloseq)
library(MCMCglmm)
library(knitr)
library(kableExtra)
library(cowplot)

source("Scripts/news_ASVs/0_Functions.R")

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

minAb=0 

# load data 
# ---------

protists_taxonomy <- read_csv("Data/Final/Protists_taxonomy_functions.csv") %>%
  mutate(Functions=ifelse(is.na(Functions),'Unknown',Functions))

pooled_OTU_table = as.matrix(read.table("Data/Final/Pooled_Protists_ASVs_counts.txt",check.names = F))

pooled_metadata = read_csv("Data/Final/Pooled_metadata.csv") %>%
  mutate(plotID=as.character(plotID),
         Alliance2 = gsub("_"," ",Alliance)) 

# Prep data: Create phyloseq object 
# ----------------------------------

MT <- as.data.frame(pooled_metadata); rownames(MT)=MT$plotID
PT <- as.matrix(protists_taxonomy); rownames(PT)=PT[,'seqs']
FCT <- as.matrix(protists_taxonomy); rownames(FCT)=FCT[,'seqs']; FCT <- FCT[,c("Functions",'PR2_Class')]

protist_PS_FCT=phyloseq(tax_table(FCT),otu_table(pooled_OTU_table,taxa_are_rows=T),sample_data(MT))

FCT_PS <-  protist_PS_FCT %>%
  tax_glom(taxrank = 'Functions' ,NArm=F) %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>% 
  mutate(Functions=factor(Functions,levels = c("Unknown","Phototrophic","Parasite","Consumer")))

Trait_Ratios_pooled <- FCT_PS %>%
  pivot_wider(id_cols = plotID, names_from = Functions,values_from = Abundance) %>%
  left_join(pooled_metadata,by="plotID") %>%
  mutate(Ratio_Parasite_Consumer=Parasite/Consumer,
         Ratio_Consumer_Parasite=1/Ratio_Parasite_Consumer,
         logRatio_Consumer_Parasite=log(Ratio_Consumer_Parasite,10),
         logRatio_Parasite_Consumer=log(Ratio_Parasite_Consumer,10),
         log_SoilWaterC = log(bulkSoilWaterContent),
         altitude = altitude)


summary(Trait_Ratios_pooled$Ratio_Consumer_Parasite,)
hist(Trait_Ratios_pooled$Ratio_Consumer_Parasite)


Data = Trait_Ratios_pooled %>% 
  subset(!is.na(Alliance))  %>%
  subset(!logRatio_Parasite_Consumer==-Inf)  %>%
  select(Parasite,Consumer,Unknown,Phototrophic,Ratio_Consumer_Parasite,logRatio_Consumer_Parasite,bio1,bio12,altitude,log_SoilWaterC,pH,CaO,Alliance,C.N,TOC)



# Basic stats 
# ---------
alt_bin = Data  %>% 
  mutate(altitude_bin1=cut_width(altitude, width=1000, boundary=0)) %>%
  group_by(altitude_bin1) %>% 
  summarise(mean_Ratio_Consumer_Parasite=mean(Ratio_Consumer_Parasite),
            sd_Ratio_Consumer_Parasite=sd(Ratio_Consumer_Parasite),
            median_Ratio_Consumer_Parasite=median(Ratio_Consumer_Parasite))

alt_bin

# Run analysis
# ------------

# Cross correlation between variables 
#pdf(file = "Redaction/V3/Figures/Cross_correlation.pdf")
chart.Correlation(select(Data, altitude,log_SoilWaterC, pH,CaO,C.N,TOC,bio1,bio12), histogram=TRUE, pch=19)
#dev.off()

# Altitude model - univariate 
baye_model_linear = MCMCglmm(fixed=logRatio_Consumer_Parasite ~  altitude , 
                             data=Data,
                             nitt=100000)

Data = Data %>% 
  mutate(predicted_MCMCGlmm_linear=predict.MCMCglmm(baye_model_linear),
         residual_MCMCGlmm_linear=logRatio_Consumer_Parasite-predicted_MCMCGlmm_linear) 

Fig4a = Data %>% 
  ggplot(aes(y=Ratio_Consumer_Parasite,x=altitude))+ geom_point() + 
  geom_smooth(method="lm",se=F,color="grey") + MyTheme +
  ylab("Consumer / Parasite (read counts ratio)")+ xlab("Elevation (meters ASL)") +
  scale_y_continuous(trans='log',breaks = c(1,2, 10, 50,100)) #+
  ggtitle(("A"))

post_prob = tibble(slope = baye_model_linear$Sol[,"altitude"])

Fig4b = post_prob %>%  
  ggplot(aes(x=slope*10^4)) + geom_histogram(binwidth = 0.1) +
  geom_density(aes(y= 0.1*..count..)) +
  MyTheme + theme(axis.line.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()) + 
  xlab("Slope coef. estimate (x10,000) ") +ylab("") +
  xlim(c(-.2,5)) + geom_vline(xintercept=0,linetype=3) # +
  ggtitle(("B"))

Fig4b


Fig4 <-  ggdraw() +
  draw_plot(Fig4a) +
  draw_plot(Fig4b, x = 0.2, y = .65, width = .3, height = .3)+
  theme(plot.margin = margin())
Fig4

ggsave(Fig4,filename = "Redaction/V4/Figures/Fig_4.pdf",width = 88,height = 88,units = "mm")


# Multivariate model 
DataM = Data %>% 
  select(logRatio_Consumer_Parasite,altitude,pH,log_SoilWaterC,bio1 , Alliance,C.N,TOC) %>% 
  na.omit()

lm_model_full = lm(logRatio_Consumer_Parasite ~  altitude+ pH +log_SoilWaterC+bio1 + Alliance+C.N+TOC, data=DataM)
lm_model_full

output_model = drop1(lm_model_full,test = "F")

out = data_frame(output_model)
out <- out %>% 
  mutate(`factor dropped`=rownames(output_model),
         `Sum of Sq`=round(`Sum of Sq`,3),
         RSS=round(RSS,1),
         AIC=round(AIC,1),
         `F value`=ifelse(`F value`<.01,"<.01",round(`F value`,2)),
        `Pr(>F)`=ifelse(`Pr(>F)`<.01,"<.01",round(`Pr(>F)`,2)))

nice_table = kable(out,booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) 
save_kable(nice_table,"Redaction/V5/Figures/Supp_Table1.pdf")


