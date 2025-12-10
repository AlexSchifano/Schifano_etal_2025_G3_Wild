#!/bin/bash

set -uo pipefail ## not set -e because line 29 can have non existing status if all sequences are the same (ie no mutations = no lines to grep)

## Bash/R script to count mut and annotate along tree provided as a table
## Table_Node.txt has to be space delimited
cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/8-COUNT_MUTATIONS
cp ../7-ARPIP/Output/Anc.fasta .

## Need to modify the Table_Node to not take the lines with root
grep -v "root" ../7-ARPIP/Output/Table_Node.txt > Table_Node.txt


## Run this script to make a GTF file for the alignment
../0-Scripts/S8.5_blast_GTF.sh

rm Count_mut_results.txt
echo "POS_Full Gene_ID Mutation Gene_Type Variant_Freq Mut_Type Strain_Name Ancestor Effect Synonymous Transition Radical Codon_Pos Site_4x_2x CDS_pos Codon" | sed "s/ /\t/g" > Count_mut_results.txt

while IFS=' ' read comp anc
do       
  anc=${anc}
  comp=${comp}
  file=count_${anc}.fa

  echo "Bash reads anc=$anc comp=$comp file=$file"

  # Need to make a file with the name of seq to grep from msa
  echo "$anc" > extract.txt
  echo "$comp" >> extract.txt
  seqkit grep -n -f extract.txt Anc.fasta -o $file
  rm extract.txt

  Rscript ../0-Scripts/R8_count_mut_tree.R $file $anc $comp

  rm $file

done < Table_Node.txt