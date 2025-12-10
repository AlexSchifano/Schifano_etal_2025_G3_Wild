setwd("~/Documents/C.elegans_sync/Final_files/01-DATA/7-ARPIP/Output/")
library(stringr)
Node_rel = read.table("Node_rel.txt", fill = T, header = F,skip = 1)
colnames(Node_rel) = c("Name","Parent","Child.1","Child.2")

write.table(Node_rel,sep="\t",quote = F,col.names = T,row.names = F,file = "Node_rel_fixed.txt")

Node_par = Node_rel[,1:2]
write.table(Node_par,sep="\t",quote = F,col.names = F,row.names = F,file = "Node_rel_fixed_2col.txt")

wild = read.table("wild.list")
wild.list = wild[,1]

Table_Node = Node_par[(Node_par$Name %in% wild.list) == F,]
Table_VC = Node_par[Node_par$Name %in% wild.list,]

# Doing the same but for the mitotypes
clustr <- read.table("../../5-CD-HIT/Clust_No_N_clustTrimmed_MSA_fasta_merged_wild_only.fasta.clstr",fill=T)
# cleaning clustr
clustr <- clustr[clustr$V1 != "+/100.00%",1:3]
clustr$V3 <- str_replace(clustr$V3,">","")
clustr <- clustr[,c(1,3)]
clustr$V3 <- str_sub(clustr$V3, end=-4)
# putting the name of the first strain as cluster name
# and the list of isotypes next to it
for(r in 1:nrow(clustr)){
  if (clustr[r,1] == ">Cluster") { clustr[r,1] <- paste(">Cluster",clustr[r+1,2],sep="_")
  for (cnt in 1:max(as.numeric(clustr$V1),na.rm = T)) {
    if(is.na(clustr[r+cnt,1]) == F) {
      if (grepl(">",clustr[r+cnt,1]) == F) { 
        if (clustr[r+cnt,1] == "0") {clustr[r,2] <- clustr[r+cnt,2]} else {clustr[r,2] <- paste(clustr[r,2],clustr[r+cnt,2],sep=", ")} 
      } else { break }
    }
    
  }
  }
}
clustr_final <- clustr[grepl(">",clustr$V1),]
clustr_final$V1 <- str_replace(clustr_final$V1,">Cluster_","")
colnames(clustr_final) <- c("Mitotype","Isotypes")

## Adding the isotypes into the Table_VC
# This makes it so we have first column being the actual isotype used and all the other columns are the isotypes within the cluster
Iso_split <- data.frame(str_split(clustr_final$Isotypes,pattern=",",simplify = T))
Iso_list <- data.frame()
for (ro in 1:nrow(Iso_split)) {
  for (co in 2:ncol(Iso_split)) {
    if (Iso_split[ro,co] != "") { 
      Iso_list[nrow(Iso_list)+1,1] <- Iso_split[ro,1]
      Iso_list[nrow(Iso_list),2] <- Iso_split[ro,co] }
  }
}
rm(Iso_split)
# I should now have in Iso_list the lines I am missing from Table_VC, I now need to add the Ancestor info
colnames(Iso_list) <- c("Mitotype","Isotype")
Iso_list$Parent <- ""
for (r in 1:nrow(Iso_list)){
  mito <- Iso_list[r,"Mitotype"]
  iso <- Iso_list[r,"Isotype"]
  Iso_list[r,"Parent"] <- Table_VC$Parent[Table_VC$Name == mito]
}
# Now I have the ancestor, I don't need the first column
Iso_list <- Iso_list[,2:3]
# I bind that previous table's two last columns with Table_VC to have the complete table with
colnames(Iso_list) <- colnames(Table_VC)
Table_VC_ALL <- rbind(Iso_list,Table_VC)
# Some white spaces have been introduces, remove them
Table_VC_ALL$Name <- str_replace(Table_VC_ALL$Name," ","")
Table_VC_ALL$Parent <- str_replace(Table_VC_ALL$Parent," ","")

write.table(Table_Node,sep="\t",quote = F,col.names = F,row.names = F,file = "Table_Node.txt")
write.table(Table_VC,sep="\t",quote = F,col.names = F,row.names = F,file = "Table_VC.txt")
write.table(Table_VC_ALL,sep="\t",quote = F,col.names = F,row.names = F,file = "Table_VC_ALL.txt")

