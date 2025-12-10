#!/bin/bash

## Script to run ARPIP with the outputes tree by IQ-TREE
## Using tree XZ1516

## But first, the tree has to be rooted

cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/7-ARPIP/Input

Tree=XZ1516_Clustered_Trimmed_MSA_fasta_merged_wild_only.fasta.treefile
cp ../../6-IQ-TREE/$Tree .

Seq=XZ1516_Clustered_Trimmed_MSA_fasta_merged_wild_only.fasta
cp ../../6-IQ-TREE/$Seq .

## The tree made by IQTREE is not rooted
## Rscript using ape to root the tree

Rscript ../../0-Scripts/R7_root_tree.R 

cp ../../0-REF/conf.txt .

if [[ -f conf.txt ]] ; then
	echo "Running ARPIP with conf.txt"
	ARPIP params=conf.txt
else
	echo "ERROR you have to make the conf.txt file to execute ARPIP"
	echo "See https://acg-team.github.io/bpp-arpip/arpip_ancestral_sequence_reconstruction_under_poisson_indel_proccess_inference_examples.html for more info"
	exit 1
fi

## Script to handle the output of ARPIP and turn it into a table of 2 col used to count mutations along the tree
grep ">" ../../6-IQ-TREE/XZ1516_Clustered_Trimmed_MSA_fasta_merged_wild_only.fasta | sed 's/>//g' > ../Output/wild.list
Rscript ../../0-Scripts/R7.5_Mod_Node_rel.R 
