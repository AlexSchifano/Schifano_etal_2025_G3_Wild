setwd("/Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/9-VARIANT_CALLING/ANNOTATION/")
library(stringr)

blast.out = read.table("blast.out")
genes=c("P","V","ND6","ND4L","W","E","12S_rRNA","S1","N","Y","ND1","ATP6","K","L1","S2","ND2","I","R","Q","F","CTB1","L2","COX3","T","ND4","COX1","C","M","D","G","COX2","H","16S_rRNA","ND3","ND5","A")
gene_type=c("tRNA","tRNA","CDS","CDS","tRNA","tRNA","rRNA","tRNA","tRNA","tRNA","CDS","CDS","tRNA","tRNA","tRNA","CDS","tRNA","tRNA","tRNA","tRNA","CDS","tRNA","CDS","tRNA","CDS","CDS","tRNA","tRNA","tRNA","tRNA","CDS","tRNA","rRNA","CDS","CDS","tRNA")
GTF = cbind("MtDNA",blast.out[,9:10],genes,gene_type)

write.table(GTF,file = "GTF_blast.gtf",quote=F,sep="\t",col.names = F,row.names = F)
