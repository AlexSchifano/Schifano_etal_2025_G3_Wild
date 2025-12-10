setwd("~/Documents/C.elegans_sync/Final_files/01-DATA/4.5_TRIM/")
library(seqinr)

## read the MSAs
MSA_2O = read.fasta("../4-ALIGN/MSA_fasta_merged_wild_2O.fasta",forceDNAtolower = F)
MSA_wild_only = read.fasta("../4-ALIGN/MSA_fasta_merged_wild_only.fasta",forceDNAtolower = F)

## read the blast outputs
blast_ATP6_2O = read.table("blast_ATP6_CDS_N2_2O.out")
blast_K_tRNA_2O = read.table("blast_K_tRNA_N2_2O.out")
blast_Last_tRNA_2O = read.table("blast_Last_tRNA_N2_2O.out")

blast_ATP6_wild = read.table("blast_ATP6_CDS_N2_wild_only.out")
blast_K_tRNA_wild = read.table("blast_K_tRNA_N2_wild_only.out")
blast_Last_tRNA_wild = read.table("blast_Last_tRNA_N2_wild_only.out")

## The pos of the seq in the N2 aligned seq are V9 = START V10 = END
# I want : END of ATP6 (without STOP), START of K tRNA and END of Last tRNA
ATP6_END_2O = as.numeric(blast_ATP6_2O[1,10])
K_tRNA_START_2O = as.numeric(blast_K_tRNA_2O[1,9])
Last_tRNA_END_2O = as.numeric(blast_Last_tRNA_2O[1,10])

ATP6_END_wild = as.numeric(blast_ATP6_wild[1,10])
K_tRNA_START_wild = as.numeric(blast_K_tRNA_wild[1,9])
Last_tRNA_END_wild = as.numeric(blast_Last_tRNA_wild[1,10])

## Make pos vector with pos to keep
CUT_2O = c(seq(from=1,to=ATP6_END_2O),seq(from=K_tRNA_START_2O,to=Last_tRNA_END_2O))
CUT_wild = c(seq(from=1,to=ATP6_END_wild),seq(from=K_tRNA_START_wild,to=Last_tRNA_END_wild))

## Cut MSAs
MSA_2O_CUT = list()
for (strain in 1:length(MSA_2O)) {
  seq_2O <- MSA_2O[[strain]]
  seq_cut_2O <- list(seq_2O[CUT_2O])
  names(seq_cut_2O) <- names(MSA_2O[strain])
  MSA_2O_CUT[strain] = seq_cut_2O
}
names(MSA_2O_CUT) = names(MSA_2O)


MSA_wild_only_CUT = list()
for (strain in 1:length(MSA_wild_only)) {
  seq_wild <- MSA_wild_only[[strain]]
  seq_cut_wild <- list(seq_wild[CUT_wild])
  names(seq_cut_wild) <- names(MSA_wild_only[strain])
  MSA_wild_only_CUT[strain] = seq_cut_wild
}
names(MSA_wild_only_CUT) = names(MSA_wild_only)

## Write final file
write.fasta(MSA_2O_CUT, names = names(MSA_2O_CUT), file.out = "Trimmed_MSA_fasta_merged_wild_2O.fasta")
write.fasta(MSA_wild_only_CUT, names = names(MSA_wild_only_CUT), file.out = "Trimmed_MSA_fasta_merged_wild_only.fasta")


