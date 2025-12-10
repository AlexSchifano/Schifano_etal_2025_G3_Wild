#!/bin/bash

## Script to blast the REF seq of 3 genes to get pos to cut MSA

cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA

## Grep one seq from the alignment to only have 1 blast output 
## In CD-HIT CB4932 is clustered with N2
seqkit grep -p CB4932 4-ALIGN/MSA_fasta_merged_wild_2O.fasta -o 4.5_TRIM/N2_2O.fasta
seqkit grep -p CB4932 4-ALIGN/MSA_fasta_merged_wild_only.fasta -o 4.5_TRIM/N2_wild_only.fasta

## replace gaps by N for blastn
sed -i '' 's/-/N/g' 4.5_TRIM/N2_2O.fasta
sed -i '' 's/-/N/g' 4.5_TRIM/N2_wild_only.fasta

## blast the ref cds (from NIH) onto the AB1 from Anc.fasta
blastn -query 0-REF/REF_ATP6_CDS_NOSTOP.fasta -outfmt 6 -subject 4.5_TRIM/N2_2O.fasta > 4.5_TRIM/blast_ATP6_CDS_N2_2O.out
blastn -query 0-REF/REF_K_tRNA.fasta -outfmt 6 -subject 4.5_TRIM/N2_2O.fasta > 4.5_TRIM/blast_K_tRNA_N2_2O.out
blastn -query 0-REF/REF_Last_tRNA.fasta -outfmt 6 -subject 4.5_TRIM/N2_2O.fasta > 4.5_TRIM/blast_Last_tRNA_N2_2O.out

blastn -query 0-REF/REF_ATP6_CDS_NOSTOP.fasta -outfmt 6 -subject 4.5_TRIM/N2_wild_only.fasta > 4.5_TRIM/blast_ATP6_CDS_N2_wild_only.out
blastn -query 0-REF/REF_K_tRNA.fasta -outfmt 6 -subject 4.5_TRIM/N2_wild_only.fasta > 4.5_TRIM/blast_K_tRNA_N2_wild_only.out
blastn -query 0-REF/REF_Last_tRNA.fasta -outfmt 6 -subject 4.5_TRIM/N2_wild_only.fasta > 4.5_TRIM/blast_Last_tRNA_N2_wild_only.out

Rscript 0-Scripts/R4_cut_MSA.R

rm 4.5_TRIM/N2_2O.fasta 4.5_TRIM/N2_wild_only.fasta
rm 4.5_TRIM/blast_*