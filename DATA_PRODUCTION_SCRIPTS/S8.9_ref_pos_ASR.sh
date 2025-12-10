#!/bin/bash

cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/8-COUNT_MUTATIONS/

LIST_ASR=/Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/8-COUNT_MUTATIONS/Count_Mut_results.txt
MS=/Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/8-COUNT_MUTATIONS/Anc.fasta
REF=/Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/0-REF/WS283_MtDNA.fasta
PATH_R=/Users/alesc399/Documents/C.elegans_sync/Final_files/02-ANALYSES/Scripts

# Prepare the file with the updated info
echo "POS_Full	Gene_ID	Mutation	Gene_Type	Variant_Freq	Mut_Type	Strain_Name	Ancestor	Effect	Synonymous	Transition	Radical	Codon_Pos	Site_4x_2x	CDS_pos	Codon	POS_REF" > /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/8-COUNT_MUTATIONS/Count_Mut_results_REF_POS.txt

while read p; do

	# Echo the line to extract info
	echo $p > temp.txt
	POS=$(awk '{print $1}' temp.txt) 
	STRAIN=$(awk '{print $9}' temp.txt)
	ANC=$(awk '{print $10}' temp.txt)

	# Echo the sequence to extract it from the MS file
	echo $ANC > extract.txt
	seqkit grep -n -f extract.txt $MS -o ${ANC}.fa

	# Grab a region arround the POS in order to blast
	seqkit subseq -r $(($POS - 10)):$((POS + 50)) ${ANC}.fa -o ${ANC}_mut.fasta

	# Blast the region against the REF
	blastn -query ${ANC}_mut.fasta -subject $REF -task blastn -outfmt 6 > blast_mut.txt

	# Add a check if the blast file is empty, if it is, expend the region searched for
	if [ -s blast_mut.txt ] ; then

		# Grab the best hit start pos
		REF_POS10=$(awk 'FNR == 1 {print $9}' blast_mut.txt)

		# Grab the query start pos
		Q_START=$(awk 'FNR == 1 {print $7}' blast_mut.txt)
		# Removing one so that the final measure is correct
		Q_START_CORREC=$(($Q_START - 1))

		# Add 10 to get the actual REF POS
		REF_POS=$(($REF_POS10 - Q_START_CORREC + 10))

	else

		# Echo a message
		echo "Initial narrow BLASTN search fail, widening the query sequence 200bp arround the mutation"

		# Grab a region arround the POS in order to blast
		seqkit subseq -r $(($POS - 200)):$((POS + 200)) ${ANC}.fa -o ${ANC}_mut.fasta

		# Blast the region against the REF
		blastn -query ${ANC}_mut.fasta -subject $REF -task blastn -outfmt 6 > blast_mut.txt

		# Grab the best hit start pos
		REF_POS10=$(awk 'FNR == 1 {print $9}' blast_mut.txt)

		# Grab the query start pos
		Q_START=$(awk 'FNR == 1 {print $7}' blast_mut.txt)
		# Removing one so that the final measure is correct
		Q_START_CORREC=$(($Q_START - 1))

		# Add 10 to get the actual REF POS
		REF_POS=$(($REF_POS10 - Q_START_CORREC + 200))

	fi

	# Add this info to the last column of the LIST_ASR
	RScript ${PATH_R}/R_pos_ref_ASR.R $POS $STRAIN $ANC $REF_POS $LIST_ASR

	rm temp.txt
	rm extract.txt
	rm ${ANC}.fa
	rm ${ANC}_mut.fasta
	rm ${ANC}.fa.seqkit.fai
	rm blast_mut.txt

done < $LIST_ASR  


