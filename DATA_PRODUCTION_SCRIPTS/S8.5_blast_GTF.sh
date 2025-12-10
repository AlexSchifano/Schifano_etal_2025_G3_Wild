#!/bin/bash

cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/8-COUNT_MUTATIONS
## Grep one seq from the alignment to only have 1 blast output 
## In CD-HIT CB4932 is clustered with N2
seqkit grep -p CB4932 Anc.fasta -o N2.fasta

## replace gaps by N for blastn
sed -i '' 's/-/N/g' N2.fasta

## blast the ref cds (from NIH) onto the AB1 from Anc.fasta
blastn -query ../0-REF/REF_RNA+CDS.fasta -outfmt 6 -subject N2.fasta > blast_CDS_RNA_N2.out

## R script to clean the output and only get the Subject pos, making a gtf
Rscript ../0-Scripts/R8.5_blast2gtf.R

rm N2.fasta blast_CDS_RNA_N2.out
