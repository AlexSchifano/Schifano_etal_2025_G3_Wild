setwd("~/Documents/C.elegans_sync/Pre_contamination/Final_files/02-ANALYSES/")

# Reading the file with the REF pos
Var=read.table("../01-DATA/9-VARIANT_CALLING/ANNOTATION/DATA_Variants_GATK_REF_POS.txt",header=T,sep="\t")

##### List of parallel mutations #####
# Quick list of candidate positions
pos_shared <- as.data.frame(table(Var$POS_REF))
# Subset at least 2 count
pos_shared <- pos_shared[pos_shared$Freq > 1,]
# Loop to get parallel mutations
Paralel_muts <- data.frame()
for (pos in pos_shared$Var1) {
  mega <- Var[Var$POS_REF == pos,]
  for(i in 1:nrow(mega)) {
    for (f in 1:nrow(mega)) {
      if (i != f) { if(sum(mega[i,1:3] == mega[f,1:3]) == 3) {
        Paralel_muts <- rbind(Paralel_muts,mega[f,])
      }}}}}

# removing natural repeats coming from the loop
Paralel_muts <- unique(Paralel_muts)

# Count for each position in Paralel_muts
Strain_shared_mut <- data.frame(table(Paralel_muts$Strain_Name))
colnames(Strain_shared_mut) <- c("Strain","NB_Shared_Mut")
write.table(Strain_shared_mut,file="List_strains_with_shared_mutations.txt",sep="\t",quote = F,row.names = F)
Strain_shared_mut <- read.table("List_strains_with_shared_mutations.txt",sep="\t",header = T)

# List of strains with suspected contamination (manual inspection)
conta_list <- c("ECA2659","JU2587","JU346","ECA1186","ECA1193","JU2519","JU2572",
              "QX1233","JU406","JU3128","NIC272","JU2825","JU2316","JU2619","JU393",
              "JU394","ECA2452","JU3127","XZ2020")

# Subset the mutations table table
Strain_Conta <- Paralel_muts[Paralel_muts$Strain %in% conta_list,]

# Subset the mutation count table
count_conta <- Strain_shared_mut[Strain_shared_mut$Strain %in% conta_list,]

# Convert Variant_Freq to numeric
Strain_Conta$Variant_Freq <- as.numeric(Strain_Conta$Variant_Freq)

# Plot with labels for mutation counts
library(ggplot2)
ggplot(Strain_Conta, aes(x = Strain_Name, y = Variant_Freq)) +
  geom_boxplot() +
  theme_minimal() +
  geom_label(data = count_conta, aes(x = Strain, y = max(Strain_Conta$Variant_Freq) + 0.1, label = NB_Shared_Mut), inherit.aes = FALSE, size = 3) +
  annotate("text",x= 10, y =  max(Strain_Conta$Variant_Freq) + 0.15, label = "Number of shared variants") +
  labs(x = "Strain", y = "Variant Frequencies")





