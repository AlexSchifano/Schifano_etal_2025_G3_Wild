#!/bin/bash
cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA
## alignement

mafft --auto --preservecase 3-MS_FASTA/fasta_merged_wild_2O.fasta > 4-ALIGN/MSA_fasta_merged_wild_2O.fasta
mafft --auto --preservecase 3-MS_FASTA/fasta_merged_wild_only.fasta > 4-ALIGN/MSA_fasta_merged_wild_only.fasta