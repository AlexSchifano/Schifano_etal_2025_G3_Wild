#!/bin/bash

# CD-HIT script, clusters sequences from fasta

MSA2O=Trimmed_MSA_fasta_merged_wild_2O.fasta
MSA=Trimmed_MSA_fasta_merged_wild_only.fasta

cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA
# cd-hit doesnt like '-' gaps, so need to replace '-' by 'N'
## MSA2O
sed 's/-/N/g' 4.5_TRIM/$MSA2O > 5-CD-HIT/No_N_${MSA2O}
cd-hit-est -i 5-CD-HIT/No_N_${MSA2O} -o 5-CD-HIT/No_N_clust_${MSA2O} -c 1

grep ">" 5-CD-HIT/No_N_clust_${MSA2O} | sed 's/>//g' > 5-CD-HIT/clust_2O.list
rm 5-CD-HIT/No_N_${MSA2O} 5-CD-HIT/No_N_clust_${MSA2O}

seqkit grep -f 5-CD-HIT/clust_2O.list 4.5_TRIM/$MSA2O -o 5-CD-HIT/Clustered_${MSA2O}


## MSA wild only
sed 's/-/N/g' 4.5_TRIM/$MSA > 5-CD-HIT/No_N_${MSA}
cd-hit-est -i 5-CD-HIT/No_N_${MSA} -o 5-CD-HIT/Clust_No_N_clust${MSA} -c 1

grep ">" 5-CD-HIT/Clust_No_N_clust${MSA} | sed 's/>//g' > 5-CD-HIT/clust_wild.list
rm 5-CD-HIT/No_N_${MSA} 5-CD-HIT/Clust_No_N_clust${MSA}

seqkit grep -f 5-CD-HIT/clust_wild.list 4.5_TRIM/$MSA -o 5-CD-HIT/Clustered_${MSA}

echo "Clustered Seq MSA can be found as 5-CD-HIT/Clustered_${MSA2O} and 5-CD-HIT/Clustered_${MSA}"
## check 
#grep -c ">" 5-CD-HIT/Clustered_Trimmed_MSA_fasta_merged_wild_*

#echo "The MSA has been clustered, see output in 5-CD-HIT/"
# cd-hit : cluster peptide sequence
# cd-hit-est : cluster nucleotide sequence
