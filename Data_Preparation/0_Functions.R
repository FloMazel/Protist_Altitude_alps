
# --------------------------------
# extraction functions
# --------------------------------

#taxo_rank="PR2_Class"
#taxo_group="Apicomplexa" 
#physeq_object=protist_PS
#rar=0
#minRead=400

subset_physeqPR2_class= function(taxo_rank,taxo_group,minRead,rar,physeq_object){
  #taxo_group = 'Parasite'
  #physeq_object =  protist_PS
  #print(minRead)
  
  oldMA <- as(tax_table(physeq_object), "matrix")
  oldDF <- data.frame(oldMA)
  newDF <- subset(oldDF, PR2_Class==taxo_group)
  newMA <- as(newDF, "matrix")
  
  new_physeq_object=physeq_object
  tax_table(new_physeq_object) <- tax_table(newMA)
  
  #subset if rarefaction
  otu_table = otu_table(new_physeq_object)
  depth = tibble(plotID=colnames(otu_table),depth=apply(otu_table,2,sum))
  
  newPlot = depth %>% 
    subset(depth>rar) %>% 
    pull(plotID)
  
  
  oldDF <- as(sample_data(new_physeq_object), "data.frame")
  newDF <- subset(oldDF, 
                  readCounts_Gregarinomorphea>minRead&
                    readCounts_Oomycota>minRead&
                    !is.na(Alliance)&
                    !is.na(TOC)&
                    plotID%in%newPlot)
  
  sample_data(new_physeq_object) <- sample_data(newDF)
  return(new_physeq_object)
  
  
}


subset_physeq= function(taxo_rank,taxo_group,minRead,rar,physeq_object){
  #taxo_group = 'Parasite'
  #physeq_object =  protist_PS
  #print(minRead)
  
  oldMA <- as(tax_table(physeq_object), "matrix")
  oldDF <- data.frame(oldMA)
  newDF <- subset(oldDF, Functions==taxo_group)
  newMA <- as(newDF, "matrix")
  
  new_physeq_object=physeq_object
  tax_table(new_physeq_object) <- tax_table(newMA)

  #subset if rarefaction
  otu_table = otu_table(new_physeq_object)
  depth = tibble(plotID=colnames(otu_table),depth=apply(otu_table,2,sum))
 
  newPlot = depth %>% 
    subset(depth>rar) %>% 
    pull(plotID)
   
   oldDF <- as(sample_data(new_physeq_object), "data.frame")
   newDF <- subset(oldDF, 
                   readCounts_Consumer>minRead&
                    readCounts_Parasite>minRead&
                      readCounts_Phototrophic>minRead&
                        !is.na(Alliance)&
                          !is.na(TOC)&
                           !is.na(C.N)&
                            plotID%in%newPlot)
   
   sample_data(new_physeq_object) <- sample_data(newDF)
   return(new_physeq_object)
   
   
}

subset_meta= function(physeq_object){
  meta = as_tibble(sample_data(physeq_object))
  return(meta)}

retrieve_meta= function(taxo_rank,taxo_group,physeq_object){
  #taxo_rank= 'functions'
  #taxo_group = 'parasite'
  #physeq_object =  protist_PS
  taxo <- tax_table(physeq_object)
  new_taxo <- taxo[taxo[,taxo_rank]%in%taxo_group,]
  new_physeq_object=phyloseq(new_taxo,otu_table(physeq_object),sample_data(physeq_object))
  
  overall_depth <- apply(as.matrix(otu_table(physeq_object)),2,sum)
  
  return(as_tibble(sample_data( new_physeq_object)))
  
}


# --------------------------------
# Diversity estimation functions
# --------------------------------

#taxo_rank="PR2_Class"
#taxo_group="Apicomplexa" 
#physeq_object=protist_PS

estimate_richness_by_taxo_ranks = function(taxo_rank,taxo_group,physeq_object){
  
  taxo <- tax_table(physeq_object)
  new_taxo <- taxo[taxo[,taxo_rank]%in%taxo_group,]
  new_physeq_object=phyloseq(new_taxo,otu_table(physeq_object),sample_data(physeq_object))
  
  overall_depth <- apply(as.matrix(otu_table(physeq_object)),2,sum)
  
  richness_protist <- estimate_richness(otu_table(new_physeq_object,taxa_are_rows=T),measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
  richness_protist$readCounts <- apply(as.matrix(otu_table(new_physeq_object)),2,sum)
  richness_protist$relativeReadCounts <- richness_protist$readCounts / overall_depth
  
  colnames(richness_protist) = paste(colnames(richness_protist),taxo_group,sep="_")
  
  return(richness_protist)
  
}

estimate_richness_by_taxo_ranks_wide_output = function(taxo_rank,taxo_group,physeq_object){
  taxo <- tax_table(physeq_object)
  new_taxo <- taxo[taxo[,taxo_rank]%in%taxo_group,]
  new_physeq_object=phyloseq(new_taxo,otu_table(physeq_object),sample_data(physeq_object))
  richness_protist <- estimate_richness(otu_table(new_physeq_object,taxa_are_rows=T),measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
  richness_protist$depth <- apply(as.matrix(otu_table(new_physeq_object)),2,sum)
  colnames(richness_protist) = paste(colnames(richness_protist),taxo_group,sep="_")
  return(richness_protist)
}

# function that estimate richness (alpha diversity metrics) given a phyloseq object and a taxonomic rankj and groups 
estimate_Beta = function(rarefaction=F,physeq_object){

  if (!rarefaction==F){physeq_object=rarefy_even_depth(physeq_object,sample.size=rarefaction)}
  
  overall_depth <- apply(as.matrix(otu_table(physeq_object)),2,sum)
  
  # transformation by sequqncing depth (normalisation)
  relative_ab <- physeq_object %>%
    transform_sample_counts(function(x) x/sum(x))
  
  beta_protist <- phyloseq::distance(physeq_object,'bray')
  beta_protist <- as.matrix(beta_protist)
  return(beta_protist)
  
}

# function that estimate richness (alpha diversity metrics) given a phyloseq object and a taxonomic rankj and groups 
estimate_Beta_by_taxo_ranks = function(taxo_rank,taxo_group,rarefaction=F,physeq_object){
  #taxo_rank= 'functions'
  #taxo_group = 'parasite'
  #physeq_object =  protist_PS
  
  taxo <- tax_table(physeq_object)
  new_taxo <- taxo[taxo[,taxo_rank]%in%taxo_group,]
  new_physeq_object=phyloseq(new_taxo,otu_table(physeq_object),sample_data(physeq_object))
  
  if (!rarefaction==F){new_physeq_object=rarefy_even_depth(new_physeq_object,sample.size=rarefaction)}
  
  overall_depth <- apply(as.matrix(otu_table(physeq_object)),2,sum)
  
  # transformation by sequqncing depth (normalisation)
  relative_ab <- new_physeq_object %>%
    transform_sample_counts(function(x) x/sum(x))
  
  beta_protist <- phyloseq::distance(new_physeq_object,'bray')
  
  return((beta_protist))
  
}


# --------------------------------
# ordiantion and tests 
# --------------------------------

metaMDS_wrapperOLD = function(beta,meta,k){
  beta = as.matrix(beta)
  beta=as.dist(beta[meta$plotID,meta$plotID]) 
  
  ana = metaMDS(beta,k)
  tmp=ana$points
  ordination <- as_tibble(tmp) %>%
    mutate(plotID=rownames(tmp))
}


extractMDS_coord = function(ana){
  tmp=ana$points
  ordination <- as_tibble(tmp) %>%
    mutate(plotID=rownames(tmp))
}

metaMDS_wrapper = function(beta,meta,k){
  beta = as.matrix(beta)
  beta=as.dist(beta[meta$plotID,meta$plotID]) 
  ana = metaMDS(beta,k)
}

metaMDS_1 = function(beta,meta,k,minRead=99){
  beta = as.matrix(beta)
  #reads=meta$readCounts_Phototrophic>minRead&meta$readCounts_Consumer>minRead&meta$readCounts_Parasite>minRead
  #to_keep=!is.na(meta$Alliance)&apply(apply(beta,1,is.na),1,sum)==0&reads # no NA and meta ok
  meta = meta %>%
    #subset(readCounts_Phototrophic>minRead&readCounts_Consumer>minRead&readCounts_Parasite>minRead) %>%
  subset(readCounts_Consumer>minRead&readCounts_Parasite>minRead) %>%
  subset(!is.na(Alliance)& !is.na(TOC)& plotID %in% colnames(beta)) 
  
  beta=as.dist(beta[meta$plotID,meta$plotID]) 
  tmp=metaMDS(beta,k)$points
  ordination <- as_tibble(tmp) %>%
    mutate(plotID=rownames(tmp))
}


permanova_wrapper = function(beta,meta,minRead=100)
{
  beta = as.matrix(beta)
  beta1=as.dist(beta[meta$plotID,meta$plotID]) 
  
  model = adonis2(beta1 ~ altitude + log(bulkSoilWaterContent) + Alliance + pH +  TOC + C.N, by="margin",data=meta)
  results = as_tibble(model)
  results$predictor = rownames(model)
  
  return(results)
  
}

permanova_soil = function(beta,meta,minRead=100)
{
  beta = as.matrix(beta)
  
  meta1 = meta %>%
    subset(readCounts_Phototrophic>minRead&readCounts_Consumer>minRead&readCounts_Parasite>minRead) %>%
    #subset(readCounts_Consumer>minRead&readCounts_Parasite>minRead) %>%
    subset(!is.na(Alliance)& !is.na(TOC)& plotID %in% colnames(beta)) %>% 
    select(plotID,altitude,bulkSoilWaterContent,Alliance,pH,TOC)
  
  beta1=as.dist(beta[meta1$plotID,meta1$plotID]) 
  
  model = adonis2(beta1 ~ altitude + log(bulkSoilWaterContent) + Alliance + pH +  TOC, by="margin",data=meta1)
  results = as_tibble(model)
  results$predictor = rownames(model)
  
  return(results)
  
}

# --------------------------------
# PLOTS
# --------------------------------

plot_altitude_ordination_with_legend = function(all,group){
  plot =all %>%
    #subset(PCR_replicate==1)
    ggplot(aes(y=MDS1,x=MDS2,color=altitude))+geom_point()+MyTheme+ggtitle(paste('NMDS',group,sep=" "))+
    scale_color_viridis_c() +
    theme(legend.position="bottom")+
  guides(color=guide_colourbar(title="Elevation (meters ASL)",title.position = "top",title.hjust=.5))
}

plot_altitude_ordination = function(all,group){
  plot = all %>%
    #subset(PCR_replicate==1)
    ggplot(aes(y=MDS1,x=MDS2,color=altitude))+geom_point()+MyTheme+
    scale_color_viridis_c() + 
    theme(legend.position = "none")+ggtitle(" ")
}


plot_alliance_ordination = function(all,group){
  plot = all %>%
    ggplot(aes(y=MDS1,x=MDS2,color=Alliance,shape=Alliance))+geom_point()+MyTheme+
    ggtitle(group)+  
    scale_shape_manual(values = alliance_display$shape) +
    scale_color_manual(values = alliance_display$color) + theme(legend.position = "none")
}

plot_soilWaterC_ordination = function(all,group){
  plot = all %>%
    ggplot(aes(y=MDS1,x=MDS2,color=log(bulkSoilWaterContent)))+geom_point()+MyTheme+ggtitle(group) + 
    theme(legend.position = "none")
}



# --------------------------------
# DNA sequqnces manipulations
# --------------------------------
extractseqnames=function(x){
  x=str_split(x,'_F_filt')[[1]][1]
  x=substring(x, 2)
  return(x)
}

unique_NA=function(x){uni = unique(x) # funciton that return unique valus expect NAs (ignored)
return(uni[!is.na(uni)])}

dualprimerHits <- function(primer1,primer2, fnFs,fnRs) {
  # Counts number of reads in which the primer is found
  nhits1 <- vcountPattern(primer1, sread(readFastq(fnFs)), fixed = FALSE)
  nhits2 <- vcountPattern(primer2, sread(readFastq(fnRs)), fixed = FALSE)
  return(sum((nhits1 > 0)& (nhits2 > 0)))
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

get_barcodes=function(x){
  first_split=strsplit(x,"-")[[1]]
  BarcodeF=first_split[1]
  BarcodeR=strsplit(first_split[2],".",fixed=TRUE)[[1]][1]
  read=strsplit(first_split[2],".",fixed=TRUE)[[1]][2]
  return(tibble(sampleID=x,BarcodeF=BarcodeF,BarcodeR=BarcodeR,read=read))
}

import_analyze_runs=function(run, path_to_runs,splitName=F){
  
  path=paste(path_to_runs,run,sep="")
  list_samples=list.files(path)
  fastq_file_characteristics=map_dfr(list_samples,import_analyze_fastq,path,run,splitName=splitName)
  return(fastq_file_characteristics)
  
}


extract_barcode_pairs=function(fastq_files_names){
  
  first_split=strsplit(fastq_files_names,"-")[[1]]
  BarcodeF=first_split[1]
  BarcodeR=strsplit(first_split[2],".",fixed=TRUE)[[1]][1]
  F_R_barcode=paste(BarcodeF,BarcodeR,sep="_")
  read=strsplit(first_split[2],".",fixed=TRUE)[[1]][2]
  
  results=tibble(
    BarcodeF=BarcodeF,
    BarcodeR=BarcodeR,
    F_R_barcode= F_R_barcode,
    read=read,
    fastq_files_name=fastq_files_names)
  
  return(results)
}

import_analyze_fastq=function(sampleID, path_to_demultiplexed,run,splitName=F){
  
  if (splitName){
    # get barcodes by formating sample names 
    first_split=strsplit(sampleID,"-")[[1]]
    BarcodeF=first_split[1]
    BarcodeR=strsplit(first_split[2],".",fixed=TRUE)[[1]][1]
    F_R_barcode=paste(BarcodeF,BarcodeR,sep="_")
    read=strsplit(first_split[2],".",fixed=TRUE)[[1]][2]
  }
  
  
  # Read fastq file and extract features 
  fastq_file=paste(path_to_demultiplexed,sampleID,sep="/")
  imported_fastq=readFastq( fastq_file)
  fastq=imported_fastq@sread@ranges@width
  
  if (splitName==FALSE){
    results=tibble(
      fastqfile_name=sampleID,
      run=run,
      minLength=min(fastq),
      maxLength=max(fastq),
      medianLength=median(fastq),
      meanLength=mean(fastq),
      read_number=length(imported_fastq@sread)
    )
  }
  else {
    results=tibble(
      #studyID=study_name,
      #sample_name=fastq_name,
      BarcodeF=BarcodeF,
      BarcodeR=BarcodeR,
      F_R_barcode= F_R_barcode,
      read=read,
      fastqfile_name=sampleID,
      run=run,
      minLength=min(fastq),
      maxLength=max(fastq),
      medianLength=median(fastq),
      meanLength=mean(fastq),
      read_number=length(imported_fastq@sread)
    )
  }
  
  return(results)
}


# --------------
# Misc functions 
# --------------

gsub1=function(x,pattern, replacement){gsub(pattern, replacement, x)}

# retrieve a subseted vector with matching pattern 
grepl_wrap=function(pattern,x){
  x[grepl(pattern,x)]
}


filter_rare=function(name, x,raretaxa){
  OTU_table=x[[name]]
  
  prevASVs=apply(OTU_table>0,1,sum)
  to_keep=prevASVs>raretaxa
  OTU_table = OTU_table[to_keep,]
  
  
  return(OTU_table)
  
}


