setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")

###### Load the script with the functions I made ######
source("~/Documents/Scripts/R_FUNCTIONS.R")

##### Compute correlation between nb of syn and nsyn muts #####
## Load mutation table
MEGA <- read.table("R_txt_files/Mutation_List.txt",header=T,sep="\t")

## Load GTF and name the columns
V170_REF_GTF <- read.table("../01-DATA/9-VARIANT_CALLING/GTF/GTF_V170.fa.gtf")
colnames(V170_REF_GTF) <- c("Chr","START","END","Gene_ID","Gene_Type")
# Need to add a column with size
V170_REF_GTF$SIZE <- V170_REF_GTF$END - V170_REF_GTF$START + 1

#### For CDS Syn ####
{
  # Getting only CDS
  MEGA_CDS <- MEGA[MEGA$Gene_Type == "CDS",]
  # Getting only Synonymous
  MEGA_CDS_Syn <- MEGA_CDS[MEGA_CDS$Synonymous == "Synonymous",]
  # Subsetting the GTF
  V170_REF_GTF_CDS <- V170_REF_GTF[V170_REF_GTF$Gene_Type == "CDS",]
  # Function to get gene count !!! Do not use mut freq here sinc it has to be computed using numner of synonymous sites
  # See script R_Count_S_NS_Sites.R for more info
  mut_freq_gene(MEGA_CDS_Syn,V170_REF_GTF_CDS)
  # Load the counts for NSyn and Syn sites per gene
  Syn_sites_genes <- read.table("R_txt_files/Syn_NSyn_Sites_count.txt",header = T,sep="\t")
  # Using it to get the correct frequencies
  Mut_freqs_MEGA_CDS_Syn$MutationFreq <- Mut_freqs_MEGA_CDS_Syn$MutationCount / Syn_sites_genes$Syn
  Mut_freqs_MEGA_CDS_Syn$SIZE <- Syn_sites_genes$Syn
  colnames(Mut_freqs_MEGA_CDS_Syn) <- c("Gene_ID","MutationCount","START","Nb_Syn_Sites","Gene_Type","MutationFreq" )
}

##### For CDS NSyn #####
{   
  # Getting only CDS
  MEGA_CDS <- MEGA[MEGA$Gene_Type == "CDS",]
  # Getting only Synonymous
  MEGA_CDS_NSyn <- MEGA_CDS[MEGA_CDS$Synonymous == "NON_Synonymous",]
  # Subsetting the GTF
  V170_REF_GTF_CDS <- V170_REF_GTF[V170_REF_GTF$Gene_Type == "CDS",]
  # Function to get gene count !!! Do not use mut freq here sinc it has to be computed using numner of synonymous sites
  # See script R_Count_S_NS_Sites.R for more info
  mut_freq_gene(MEGA_CDS_NSyn,V170_REF_GTF_CDS)
  # Load the counts for NSyn and Syn sites per gene
  Syn_sites_genes <- read.table("R_txt_files/Syn_NSyn_Sites_count.txt",header = T,sep="\t")
  # Using it to get the correct frequencies
  Mut_freqs_MEGA_CDS_NSyn$MutationFreq <- Mut_freqs_MEGA_CDS_NSyn$MutationCount / Syn_sites_genes$NSyn
  Mut_freqs_MEGA_CDS_NSyn$SIZE <- Syn_sites_genes$NSyn
  colnames(Mut_freqs_MEGA_CDS_NSyn) <- c("Gene_ID","MutationCount","START","Nb_Syn_Sites","Gene_Type","MutationFreq" )
}

#### Combine tables ####
{
  Mut_gene <- Mut_freqs_MEGA_CDS_Syn[,-5]
  colnames(Mut_gene) <- c("Gene_ID", "Syn_Count", "START", "Nb_Syn_Sites", "Syn_Freq")
  Mut_gene$NSyn_Count <- Mut_freqs_MEGA_CDS_NSyn$MutationCount
  Mut_gene$Nb_NSyn_Sites <- Mut_freqs_MEGA_CDS_NSyn$Nb_Syn_Sites
  Mut_gene$NSyn_Freq <- Mut_freqs_MEGA_CDS_NSyn$MutationFreq
}

#### Plot ####
{
  library(ggplot2)
  
  ggplot(Mut_gene, aes(x=Syn_Freq, y=NSyn_Freq)) + geom_line()
  
}

#### Test ####
{
  lm_model <- lm(Mut_gene$Syn_Freq ~ Mut_gene$NSyn_Freq)
}







