#!/bin/bash

## Script for IQ-TREE making one tree and then another using the most basal wild strain as new outgroup

cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/6-IQ-TREE

## Tree with clustered lines rooted with brenneri and tropicalis
ALN_2O=Clustered_Trimmed_MSA_fasta_merged_wild_2O.fasta
cp ../5-CD-HIT/$ALN_2O .
iqtree2 -s $ALN_2O -m MFP -o MtDNA_C.brenneri,MtDNA_C.tropicalis

## In this first tree I can see that the cluster made of XZ1516 | ECA1493 | ECA1713 | ECA2195 | ECA2191 | ECA2199
## is basal to the rest of the wild isolates, i'll use one of those strains (XZ1516) as the reference seq for this cluster 
## and as the outgroup for a tree made with only wild isolates

## Tree with wilds only, rooted on XZ1516 ==> remove other lines from the same cluster
ALN_W=Clustered_Trimmed_MSA_fasta_merged_wild_only.fasta
cp ../5-CD-HIT/$ALN_W .
## I need to remove the lines in the same cluster as XZ1516 = ECA1493 | ECA1713 | ECA2195 | ECA2191 | ECA2199
echo "ECA1493" > sub.txt
echo "ECA1713" >> sub.txt
echo "ECA2195" >> sub.txt
echo "ECA2191" >> sub.txt
echo "ECA2199" >> sub.txt
seqkit grep -v -f sub.txt $ALN_W -o XZ1516_${ALN_W}
rm sub.txt
## Make tree rooted on XZ1516
iqtree2 -s XZ1516_${ALN_W} -m MFP -o XZ1516