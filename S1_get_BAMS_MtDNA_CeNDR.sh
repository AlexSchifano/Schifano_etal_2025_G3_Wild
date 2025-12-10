#!/bin/bash

## Script for full procedure of extracting MtDNA data from CeNDR bamfils, indexing, aligning, triming and making the final tree

# PATH on office machine
cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA

## !!! The following isotypes are not available as diplayed name : CB4851 CB4853 CB4855 CB4857 CB4858 PB306 check CaeNDR to download
## EDIT: Try and replace by ECA243 ECA246 ECA248 ECA250 ECA251 ECA259 (matching order of previous line)

# list of isotype references in isotype.list, see Technical stuff Word file for details
# isotype.list does not contain the following lines:
# ECA2659 ; JU2587 ; JU346 ; ECA1186 ; ECA1193 ; JU2519 ; JU2572 ; QX1233 ; JU406 ; JU3128 ; NIC272 ; JU2825 ; JU2316 ; JU2619 ; JU393 ; JU394 ; ECA2452 ; JU3127 ; XZ2020
# These lines are suspected to have mtDNA contamination from different C. elegans lines due to several variants showing very similar frequencies, which we deem to be 
# unlikely and we choose to remove those lines 
for NAME in $(cat isotype.list) ; do

	BAM=${NAME}.bam

	## 1) Download from CeNDR with the script on website https://caendr.org/data/data-release/c-elegans/latest
	# It dowloads full genome bam, with nuclear for ALL ~1700 strains !!!!!
	# Or use:

	## Check if the strain bam has already been downloaded or not
	if [[ -s "1-MtDNA_BAM/${NAME}_MtDNA.bam" ]] ; then
		echo "${NAME}_MtDNA.bam already exists"
	else 
		wget -O "$BAM" "https://storage.googleapis.com/caendr-site-private-bucket/bam/c_elegans/$BAM"

		## 2) Index the bam files
		samtools index $BAM

		## 3) extract MtDNA data
		samtools view -b $BAM MtDNA > 1-MtDNA_BAM/${NAME}_MtDNA.bam

		## remove the full bam file (too big)
		rm $BAM ${BAM}.bai 

	fi

		## get consensus from the bam | Using -m simple to not comput Ns if they are in small numbers 
		# with -c set at 0.51 to have the consensus be the nuc at least present in 51% of reads

		samtools consensus -c 0.50 -m simple -f fasta 1-MtDNA_BAM/${NAME}_MtDNA.bam  -o 2-STRAIN_FASTA/${NAME}_MtDNA.fasta
		
done

	## Align
	## Cut the bad regions (interegnic between ATP6 and K tRNA)
	## Make tree et ectera


## Reference genome WS283
#wget -O "20220216_c_elegans_WS283.genome.fa" "https://storage.googleapis.com/caendr-site-public-bucket/dataset_release/c_elegans/20220216/20220216_c_elegans_WS283.genome.fa"
#mv 20220216_c_elegans_WS283.genome.fa Elegans_WS283_Full_refgenome.fa
#mv Elegans_WS283_Full_refgenome.fa REF/

## Need to provide a bed file stating which regions to cut for samtools ampliconclip:
#echo "\nMtDNA 3233 3244\nMtDNA 13328 14000" | sed "s/ /\t/g" > elegans_cut_interg_Dloop.bed
#BED=elegans_cut_interg_Dloop.bed



