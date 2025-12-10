setwd("~/Documents/C.elegans_sync/Final_files/01-DATA/0-REF/")
library(stringr)
data <- read.csv("20220216_CeaNDR_strain_data.csv")
# removes the strains added in the last release (after I started my work)
data_old <- data[data$release != 20231213,]

# Putting the strains for each isotype
table <- data.frame(isotype = unique(data_old$isotype))
row.names(table) <- table$isotype
table$strains <- ""
for (f in table$isotype) {
  subset <- data_old[data_old$isotype == f,]
  table[f,"strains"] <- paste(subset$strain,collapse =", ")
}

# Doing the same but for the mitotypes
clustr <- read.table("../5-CD-HIT/Clust_No_N_clustTrimmed_MSA_fasta_merged_wild_only.fasta.clstr",fill=T)
# cleaning clustr
clustr <- clustr[clustr$V1 != "+/100.00%",1:3]
clustr$V3 <- str_replace(clustr$V3,">","")
clustr <- clustr[,c(1,3)]
# putting the name of the first strain as cluster name
# and the list of isotypes next to it
for(r in 1:nrow(clustr)){
  if (clustr[r,1] == ">Cluster") { clustr[r,1] <- paste(">Cluster",clustr[r+1,2],sep="_")
    for (cnt in 1:max(as.numeric(clustr$V1),na.rm = T)) {
      if (grepl(">",clustr[r+cnt,1]) == F) { 
        if (clustr[r+cnt,1] == "0") {clustr[r,2] <- str_sub(clustr[r+cnt,2],end=-4)} else {clustr[r,2] <- str_sub(paste(clustr[r,2],clustr[r+cnt,2],sep=", "),end=-4)} 
        } else { break }
    }
    }
}
clustr_final <- clustr[grepl(">",clustr$V1),]
clustr_final$V1 <- str_replace(clustr_final$V1,">Cluster_","")
clustr_final$V1 <- str_sub(clustr_final$V1, end=-4)
colnames(clustr_final) <- c("Mitotype","Isotypes")

## Some of the isotypes in "table" are not in the study, I look at my isotype.list
isotype.list <- read.table("../isotype.list", col.names = "Isotypes")
# I put the strains for each isotype
isotype.list$Strains <- ""
for (iso in isotype.list$Isotypes) {
  isotype.list$Strains <- table$strains[table$isotype == iso]
}
## !!! The isotypes have changed, example EC243 used to be an isotypes but is now a strain in another isotype 
## I only provide the list of mitotypes and matching isotypes I have

# Table is now done
# writing it
library(openxlsx)
write.xlsx(clustr_final,file = "../0-REF/List_Mitotype_Isotype_Strains.xlsx")






