#!/bin/bash

cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA

## Merge fasta from consensus into one fasta for alignment

ls 2-STRAIN_FASTA/ > fasta.list
rm 3-MS_FASTA/fasta_merged_wild_only.fasta

## Get fasta seq for C.briggsae and C.tropicalis, downloaded from NCBI
HEAD_B=$(head -1 0-REF/Brenneri_NC035244_Mt_REF.fasta)
HEAD_T=$(head -1 0-REF/Tropicalis_NC025756_Mt_REF.fasta)
grep "" 0-REF/Brenneri_NC035244_Mt_REF.fasta | sed "s/${HEAD_B}/>MtDNA:C.brenneri/g" > 3-MS_FASTA/fasta_merged_wild_2O.fasta
grep "" 0-REF/Tropicalis_NC025756_Mt_REF.fasta | sed "s/${HEAD_T}/>MtDNA:C.tropicalis/g" >> 3-MS_FASTA/fasta_merged_wild_2O.fasta

for fasta in $(cat fasta.list) ; do
	strain=$(basename $fasta _MtDNA.fasta)
	sed "s/>MtDNA/>${strain}/g" 2-STRAIN_FASTA/$fasta >> 3-MS_FASTA/fasta_merged_wild_2O.fasta
	sed "s/>MtDNA/>${strain}/g" 2-STRAIN_FASTA/$fasta >> 3-MS_FASTA/fasta_merged_wild_only.fasta
done

echo "Fasta are merged in 2 files, 3-MS_FASTA/fasta_merged_wild_2O.fasta has all wild isolate + 2 Outgroups and has all wild isolates 3-MS_FASTA/fasta_merged_wild_only.fasta"



## Still need to align and trim repeated regions

