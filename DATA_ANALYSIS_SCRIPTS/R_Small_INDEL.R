setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")

###### Load the script with the functions I made ######
source("~/Documents/Scripts/R_FUNCTIONS.R")

## Load mutation table
MEGA <- read.table("R_txt_files/Mutation_List.txt",header=T,sep="\t")

## Subset indels
MEGA_INDEL <- MEGA[MEGA$Mut_Type != "SNP",]
library(openxlsx)
write.xlsx(MEGA_INDEL,file="R_txt_files/Indel_List.xlsx")

## Count indels from ASR (they are labeled DEL or INS whereas INDEL is the label for the variants)
nrow(MEGA_INDEL[MEGA_INDEL$Mut_Type != "INDEL",])
nrow(MEGA_INDEL[MEGA_INDEL$Mut_Type == "INDEL",])

###### Computing density per gene type #####
{
  ## Distribution of indels
  table(MEGA_INDEL$Gene_Type)
  CDS_indel_count <- nrow(MEGA_INDEL[MEGA_INDEL$Gene_Type == "CDS",])
  tRNA_indel_count <- nrow(MEGA_INDEL[MEGA_INDEL$Gene_Type == "tRNA",])
  NC_indel_count <- nrow(MEGA_INDEL[MEGA_INDEL$Gene_Type == "NON_CODING",])
  
  ## GTF
  V170_REF_GTF = read.table("../01-DATA/9-VARIANT_CALLING/GTF/GTF_V170.fa.gtf")
  V170_REF_GTF <- subset(V170_REF_GTF, select = -V1)
  V170_REF_GTF[,5] <- V170_REF_GTF[,2] - V170_REF_GTF[,1] +1
  colnames(V170_REF_GTF) = c("START","END","Gene_ID","Gene_Type","SIZE")
  
  ## Nb bp CDS and tRNA
  Size_CDS <- sum(V170_REF_GTF$SIZE[V170_REF_GTF$Gene_Type == "CDS"])
  Size_tRNA <- sum(V170_REF_GTF$SIZE[V170_REF_GTF$Gene_Type == "tRNA"])
  
  ## Loop to count non-coding pos
  coding_pos=c()
  for (f in 1:nrow(V170_REF_GTF)) {
    coding_pos <- c(coding_pos,seq(from=V170_REF_GTF[f,"START"],to=V170_REF_GTF[f,"END"]))
  }
  all_pos = seq(from=1,to=max(V170_REF_GTF$END))
  non_coding_pos = which(all_pos %in% coding_pos == F)
  length(non_coding_pos)
  Size_NC <- length(coding_pos)
  
  ## Computing INDEL proportions based on bp
  # CDS
  CDS_freq <- CDS_indel_count / Size_CDS
  tRNA_freq <- tRNA_indel_count / Size_tRNA
  NC_freq <- NC_indel_count / Size_NC
  
  CDS_prop <- CDS_freq / sum(CDS_freq,tRNA_freq,NC_freq)
  tRNA_prop <- tRNA_freq / sum(CDS_freq,tRNA_freq,NC_freq)
  NC_prop <- NC_freq / sum(CDS_freq,tRNA_freq,NC_freq)
  
}

##### Looking for hotspots #####
{
  table(MEGA_INDEL$POS_REF)
}

##### Waneka data #####
{
  Waneka_indel <- read.xlsx("Previous_results/Waneka/File_S2.xlsx")
  
  ## GTF
  GTF_WS283 <- read.table("../01-DATA/0-REF/WS283_MtDNA.gtf", col.names = c("Chr","Start","End","Gene_ID","Gene_Type"))
  row.names(GTF_WS283) <- GTF_WS283$Gene_ID
  
  ## Annotate gene type in Waneka_indel
  Waneka_indel$gene_type <- ""
  Waneka_indel$gene_ID <- ""
  for (row in 1:nrow(Waneka_indel)) {
    pos <- Waneka_indel[row,"mtDNA.position"]
    for (gene in GTF_WS283$Gene_ID){
      range <- seq(from=GTF_WS283[gene,"Start"], to=GTF_WS283[gene,"End"])
      if (pos %in% range) { 
        Waneka_indel[row,"gene_type"] <- GTF_WS283[gene,"Gene_Type"]
        Waneka_indel[row,"gene_ID"] <- GTF_WS283[gene,"Gene_ID"] }}}
  
  ## Indel per type of gene
  CDS_nb <- nrow(Waneka_indel[Waneka_indel$gene_type == "CDS",])
  tRNA_nb <- nrow(Waneka_indel[Waneka_indel$gene_type == "tRNA",])
  rRNA_nb <- nrow(Waneka_indel[Waneka_indel$gene_type == "rRNA",])
  NC_nb <- nrow(Waneka_indel[Waneka_indel$gene_type == "",])
  
  ## Insertions vs Deletions
  table(Waneka_indel$indel.type[Waneka_indel$gene_type == "CDS"])
  table(Waneka_indel$indel.type[Waneka_indel$gene_type == "tRNA"])
  table(Waneka_indel$indel.type[Waneka_indel$gene_type == "rRNA"])
  table(Waneka_indel$indel.type[Waneka_indel$gene_type == ""])
  
  ## INDEL proportion
  CDS_f <- CDS_nb / 10296
  tRNA_f <- tRNA_nb / 1239
  rRNA_f <- rRNA_nb / 1650
  NC_f <- NC_nb / 138
  NC_f2 <- 1/138

  CDS_prop <- CDS_f / sum(CDS_f+tRNA_f+rRNA_f+NC_f2)
  tRNA_prop <- tRNA_f / sum(CDS_f+tRNA_f+rRNA_f+NC_f2)
  rRNA_prop <- rRNA_f / sum(CDS_f+tRNA_f+rRNA_f+NC_f2)
  NC_prop <- NC_f2 / sum(CDS_f+tRNA_f+rRNA_f+NC_f2)
}












