#####################################################################################################################
#     Readme file for the scripts associated with the analysis found in the paper “Soil protist function varies with elevation in the Swiss Alps”    #
#####################################################################################################################

# ——————————
# Overview of the scripts 
# ——————————

 The folder contains two subfolders: one (“Data_Preparation”) contains scripts dedicated to DNA data production and the other (“Figures_Statistics”) contains scripts dedicated to the statistical analysis and the production of figures. The statistical analysis and the production of figures can be carried out by using the prepared data (ASV table, ASV taxonomy and metadata) provided as additional files with the paper. 

# —————————
# Details of the scripts  
# —————————

1. Data_Preparation
# —————————

0_Functions.R contains the functions developed for the purpose of this analysis 
2_Demultiplexing.R contains the code that demultiplex the initial fastq file (from the sequencing facility) and remove barcode and primers from the reads 
3_Subset_Rename_Demultiplexed_fastq_files.R contains the code to reorganise samples, filter and trim and do some basic quality checks
4_ASV_inference_Merging.R contains the code to infer ASVs (dada2 pipeline)         
5_Combine_Sequencing_replicates_in_one_OTU_table.R contains the code to combine PCR replicates by site
6_Taxonomic_assignments.R contains the code to assign taxonomy to ASVs                   
7_ASVs_filtering_functional_annotations.R contains the code to filter and functional assign ASVs
8_Pooling_Sequencing_reps_and_SampleChoice.R contains the code to pool the PCR replicates and choose samples 


2. Data analysis and figure production 
# ————————————————

Each script refers to a specific figure of the paper and its associated statistical test. 
This part can be carried out independently from the data preparation part by using the prepared data (ASV table, taxonomy and metadata) provided as Additional files with the paper. 

Fig_1.R
Fig_2.R         
Fig_3.R
Fig_4.R
Fig_5.R
Fig_6.R    
Supp_Figure_1_2.R


