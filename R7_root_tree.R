setwd("/Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/7-ARPIP/Input/")

#install.packages("ape")
library(ape)

tree = read.tree("XZ1516_Clustered_Trimmed_MSA_fasta_merged_wild_only.fasta.treefile")

rooted = root(tree, outgroup = "XZ1516", resolve.root = TRUE)

# is.rooted(rooted)
# plot(rooted)

write.tree(rooted, file = "Rooted_XZ1516_tree.newick", digits = 10)
