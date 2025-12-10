setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")

library(seqinr)

## read the fasta with all the CDS from the REF
REF_CDS = read.fasta("../01-DATA/0-REF/REF_CDS.fasta",forceDNAtolower = F)

###### Load the script with the functions I made ######
source("~/Documents/Scripts/R_FUNCTIONS.R")

## Use the degeneracy function
degeneracy_tag_CDS(REF_CDS)
REF_CDS_Degen <- Seq_table_list

## Calculing Nsyn and Syn sites
Syn_sites_genes = data.frame(Gene=names(REF_CDS_Degen),NSyn=0,Syn=0,NSS=0)
rownames(Syn_sites_genes) = names(REF_CDS_Degen)
# Counting 2x and 3x sites as 1/3 Syn
degen_counts_all = data.frame(matrix(ncol=2,nrow=0))
colnames(degen_counts_all) <- c("Var1","Freq")
for (gene in names(REF_CDS_Degen)) {
  gene_table <- REF_CDS_Degen[[gene]]
  degen_counts <- table(gene_table$Degeneracy)
  degen_counts_1 <- data.frame(table(gene_table$Degeneracy))
  degen_counts_all <- rbind(degen_counts_all,degen_counts_1)
  if ("Threefold" %in% names(degen_counts)) {
    NSyn_count <- degen_counts[["Zerofold"]] + 2*(degen_counts[["Twofold"]]+degen_counts[["Threefold"]]) / 3
    Syn_count <- degen_counts[["Fourfold"]] + (degen_counts[["Twofold"]]+degen_counts[["Threefold"]]) / 3
  } else {
    NSyn_count <- degen_counts[["Zerofold"]] + 2*degen_counts[["Twofold"]] / 3
    Syn_count <- degen_counts[["Fourfold"]] + degen_counts[["Twofold"]] / 3
  }
  
  ratio <- NSyn_count/Syn_count
  Syn_sites_genes[gene,"Gene"] <- gene
  Syn_sites_genes[gene,"NSyn"] <- NSyn_count
  Syn_sites_genes[gene,"Syn"] <- Syn_count
  Syn_sites_genes[gene,"NSS"] <- ratio
}

View(Syn_sites_genes)

## Writing into a file
write.table(Syn_sites_genes, file = "R_txt_files/Syn_NSyn_Sites_count.txt",sep="\t",quote = F,row.names = F)
