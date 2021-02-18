# Plot taxo/function relationship!

rm(list=ls())
library(tidyverse)
library(networkD3)
library(webshot)

MyTheme=theme_classic() +theme(plot.title.position = "plot",plot.title = element_text(size=16),axis.text=element_text(size=14),axis.text.x =element_text(size=14,angle=45,hjust=1),axis.title=element_text(size=14),legend.text=element_text(size=12),panel.border = element_blank(),axis.ticks = element_line(size = 1.5),axis.ticks.length = unit(.4, "cm"))

# loading 
pooled_OTU_table = as.matrix(read.table("Data/Final/Pooled_Protists_ASVs_counts.txt",check.names = F))
ASV_stats = tibble(seqs=rownames(pooled_OTU_table),readCounts_pooled=apply(pooled_OTU_table,1,sum),prevalence_pooled=apply(pooled_OTU_table>0,1,sum))

updated_Protist_taxonomy <- read.csv("Data/Final/Protists_taxonomy_functions.csv",stringsAsFactors=FALSE) %>% 
  left_join(ASV_stats,by="seqs")

# Basic count stats 
dim(pooled_OTU_table)

updated_Protist_taxonomy %>% 
  group_by(Functions) %>% 
  summarise(nASV=n(),
            ReadCounts =sum(readCounts_pooled)) %>% 
  mutate(rel_nASV=nASV/sum(nASV),
         relReadCounts =ReadCounts/sum(ReadCounts))

sum(updated_Protist_taxonomy$readCounts_pooled)

# sankey diagram 

#format data
data = updated_Protist_taxonomy %>% 
  group_by(PR2_Class,Functions) %>%
  summarise(PR2_Division=unique(PR2_Division),
            PR2_Supergroup=unique(PR2_Supergroup),
            Functions=unique(Functions),
            nASVS=n(),
            val=1,
            TotalCOunts=sum(readCounts_pooled)) 

# Create link dataframe
Links = data %>% subset(TotalCOunts>1000) %>%
  mutate(Functions=ifelse(is.na(Functions),'Unknown_function',Functions),
         Taxonomy_Class=ifelse(grepl('Unknown',PR2_Class),'Unknown_class',PR2_Class), 
         Taxonomy_Division=ifelse(grepl('Unknown',PR2_Division),'Unknown_Division',PR2_Division),
         Taxonomy_Supergroup=ifelse(grepl('Unknown',PR2_Supergroup),'Unknown_Supergroup',PR2_Supergroup)) %>%
 # arrange(PR2_Supergroup) # reorder Classes by Division 
 arrange(PR2_Supergroup,PR2_Division) # reorder Classes by Division 
  
# Create node dataframe

fct = c("Phototrophic","Parasite","Consumer","Unknown_function")
fct = c("Parasite","Consumer","Phototrophic","Unknown_function")

noms = c(unique(Links$PR2_Class),fct)
Nodes_ini = data_frame(node=c(0:(length(noms)-1)),
                   names=noms)

Links = Links %>% 
  left_join(Nodes_ini,by=c('PR2_Class'='names')) %>% 
  left_join(Nodes_ini,by=c('Functions'='names')) 

# Add a 'group' column to the nodes data frame:
Nodes <- Nodes_ini %>% 
  left_join(data,by=c("names"="PR2_Class")) %>% 
  mutate(group=ifelse(names=='Phototrophic','Phototrophic',
                      ifelse(names=='Parasite','Parasite',
                             ifelse(names=='Consumer','Consumer',Functions)))) %>% 
  mutate(group=ifelse(is.na(group),"Unknown",group)) %>% 
  subset(!((names=="Chrysophyceae")&Functions%in%c('Phototrophic','Consumer'))) %>% # only keep one node ofr this group, with no functions (as it is a mix of different functions)
  mutate(names=ifelse(TotalCOunts<15000,"",names))
  
# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain(["Phototrophic", "Parasite","Consumer","Unknown","Unknown_function"]) .range(["#FFD92F", "#1F78B4","#A6761D","#BEBEBE","#BEBEBE"])'


# Create graph
Fig1a = sankeyNetwork(Links = Links, 
                      Nodes = Nodes,
                      Source = 'node.x', 
                      Target = 'node.y', 
                      NodeID = 'names',
                      Value = 'TotalCOunts',
                      colourScale=my_color, 
                      NodeGroup="group",
                      LinkGroup = "Functions",
                      iterations = 0,
                      fontSize=c(38),
                      nodePadding=5,
                      height = 2600,
                      width = 1500,
                      nodeWidth = 50,
                      fontFamily = "Arial")
Fig1a
pathNet = "/Users/fmazel/Desktop/Work/Cloud_analysis/Alpine_soil/Redaction/V3/Figures/Fig1_pooled/Fig1.html"
saveNetwork(Fig1a, pathNet, selfcontained = TRUE)

#NOW modification by hand to re-adjust position of functions and save as "Fig1_ajusted"
pathNet1 = "/Users/fmazel/Desktop/Work/Cloud_analysis/Alpine_soil/Redaction/V3/Figures/Fig1_pooled/Fig1_ajusted.htm"

webshot(pathNet1  , "Redaction/V3/Figures/Fig1_pooled/Fig1_raw.pdf")

