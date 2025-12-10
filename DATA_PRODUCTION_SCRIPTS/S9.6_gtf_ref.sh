#!/bin/bash

## Script to make a gtf for each ref (node) used for VC
cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/9-VARIANT_CALLING/ANNOTATION

REF_LOC=/Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/9-VARIANT_CALLING/REF

ls $REF_LOC | grep -v ".fai" | grep -v ".pac" | grep -v ".sa" | grep -v ".bwt" | grep -v ".amb" | grep -v ".ann" | grep -v ".dict" > Anc.list


if [[ -g ../GTF ]] ; then
	mkdir ../GTF
fi


for anc in $(cat Anc.list) ; do
	blastn -word_size 20 -query ../../0-REF/REF_RNA+CDS.fasta -outfmt 6 -subject ${REF_LOC}/$anc > blast.out

	Rscript ../../0-Scripts/R9.6_blast2gtf.R
	mv GTF_blast.gtf ../GTF/GTF_${anc}.gtf
	rm blast.out ;
done
