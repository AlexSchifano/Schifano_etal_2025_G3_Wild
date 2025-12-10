#!/bin/bash

cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA
S_PATH=/Users/alesc399/Documents/C.elegans_sync/Final_files/DATA/0-Scripts
chmod a+x 0-Scripts/*

## STEP #1 and #2, download bam/extract MtDNA reads and consensus
if [[ -d "1-MtDNA_BAM/" ]] && [[ -d "2-STRAIN_FASTA/" ]] ; then
	#1 Download files, extract MtDNA and get consensus seq
	echo "Executing Script ${S_PATH}/S1_get_BAMS_MtDNA_CeNDR.sh"
	${S_PATH}/S1_get_BAMS_MtDNA_CeNDR.sh
else
	echo "Creating Directories 1-MtDNA_BAM AND 2-STRAIN_FASTA"
	mkdir 1-MtDNA_BAM 2-STRAIN_FASTA
	#1 Download files, extract MtDNA and get consensus seq
	echo "Executing Script ${S_PATH}/S1_get_BAMS_MtDNA_CeNDR.sh"
	${S_PATH}/S1_get_BAMS_MtDNA_CeNDR.sh
fi 		

echo "STEP #1 -> #2 DONE"
echo "MtDNA bam files downloaded in 1-MtDNA_BAM/ AND consensus fasta files in 2-STRAIN_FASTA/"

## STEP #3, merge fasta in MS file
if [[ -d "3-MS_FASTA/" ]] ; then
	#2 Merge fasta
	echo "Executing Script ${S_PATH}/S2_merge_fasta.sh"
	${S_PATH}/S2_merge_fasta.sh
else 
	echo "Creating Directory 3-MS_FASTA"
	mkdir 3-MS_FASTA
	#2 Merge fasta
	echo "Executing Script ${S_PATH}/S2_merge_fasta.sh"
	${S_PATH}/S2_merge_fasta.sh
fi

echo "STEP #3 DONE"
echo "Multi Sequence FASTA files in 3-MS_FASTA, NOT ALIGNED"

## STEP #4, ALIGN MS files
if [[ -d "4-ALIGN" ]] ; then
	#3 Align
	echo "Executing Script ${S_PATH}/S3_align_mafft.sh"
	${S_PATH}/S3_align_mafft.sh
else
	echo "Creating Directory 4-ALIGN"
	mkdir 4-ALIGN
	#3 Align
	echo "Executing Script ${S_PATH}/S3_align_mafft.sh"
	${S_PATH}/S3_align_mafft.sh
fi

echo "STEP #4 DONE"
echo "MS files are Aligned in 4-ALIGN/MSA_fasta..."

## STEP #4.5, TRIMM MSA
if [[ -d "4.5_TRIM" ]] ; then
	echo "Executing Script ${S_PATH}/S4_blast_REF_trim.sh"
	${S_PATH}/S4_blast_REF_trim.sh
else
	echo "Creating Directory 4.5_TRIM"
	mkdir 4.5_TRIM
	echo "Executing Script ${S_PATH}/S4_blast_REF_trim.sh"
	${S_PATH}/S4_blast_REF_trim.sh
fi

echo "STEP #4.5 DONE"
echo "MSA files are trimmed in 4.5_TRIM/Trimmed_MSA_..."

## STEP #5 CD-HIT
if [[ -d "5-CD-HIT" ]] ; then
	echo "Executing Script S5_CD-HIT.sh"
	${S_PATH}/S5_CD-HIT.sh
else
	echo "Creating Directory 5-CD-HIT"
	mkdir 5-CD-HIT
	echo "Executing Script S5_CD-HIT.sh"
	${S_PATH}/S5_CD-HIT.sh
fi

echo "STEP #5 DONE"
echo "Sequences have been clustered in 5-CD-HIT as Clustered..."

## STEP #6 IQ-TREE
if [[ -d "6-IQ-TREE" ]] ; then
	echo "Executing Script S6_IQ-TREE.sh"
	${S_PATH}/S6_IQ-TREE.sh
else
	echo "Creating Directory 6-IQ-TREE"
	mkdir 6-IQ-TREE
	echo "Executing Script S6_IQ-TREE.sh"
	${S_PATH}/S6_IQ-TREE.sh
fi

echo "STEP #6 DONE"
echo "Both trees can be found in 6-IQ-TREE/"

## STEP #7 ARPIP
if [[ -d "7-ARPIP" ]] && [[ -d "7-ARPIP/Input" ]] && [[ -d "7-ARPIP/Output" ]]; then
	echo "Executing Script S7_ARPIP.sh"
	${S_PATH}/S7_ARPIP.sh
else
	echo "Creating Directory 7-ARPIP"
	mkdir 7-ARPIP 7-ARPIP/Input 7-ARPIP/Output
	echo "Executing Script S7_ARPIP.sh"
	${S_PATH}/S7_ARPIP.sh
fi

echo "STEP #7 DONE"
echo "ARPIP Outputs can be found in 7-ARPIP/Output"
echo "Anc.fasta has the alignment from all wild isolates and Node sequences named V.."
echo "Table_Node.txt and Table_VC.txt are also used for next steps"

## STEP #8 COUNT TREE MUTATIONS
if [[ -g "8-COUNT_MUTATIONS" ]] ; then
	echo "Executing Script S8_count_mut_tree.sh"
	${S_PATH}/S8_count_mut_tree.sh
else
	echo "Creating Directory 8-COUNT_MUTATIONS"
	mkdir 8-COUNT_MUTATIONS
	echo "Executing Script S8_count_mut_tree.sh"
	${S_PATH}/S8_count_mut_tree.sh
fi

echo "STEP #8 DONE"
echo "All results from counting mutations in the tree can be found in Count_mut_results.txt, best read in R"


## STEP #9 VARIANT CALLING + ANNOTATION
# Before run the script, the ARPIP/Output/Table_VC has to be edited.
# Table_VC.txt contains relationships between each leaf of the tree and its closest node parent.
# BUT it only includes the strains after they have been clustered by 100% match of the consensus sequence
# For the variant calling we want all 550 strains, we have to add relationships between the wild strains of a same cluster and there parents
# That way we do the Variant Calling on all wild isolates (except the 5 lines used for rooting the tree)
# This table can be found in 0-Scripts/Table_VC_All.txt

# For this step, run S9_Variant_Calling on UPPMAX

## STEP #9.5 ANNOTATION
# Once Variant Calling is done download the 9-Variant_CALLING/ Folder
# Annotation is done in 3 steps
# 9.5 = group all vcf results into one table using S9.5_vcf_results.sh and R9.5_vcf_results.R
# 9.6 = make a gtf file for each of the REF sequence (Node) using S9.6_gtf_ref.sh and R9.6_blast2gtf
# 9.9 = Final step, annotation of variants and outputs a nice table
if [[ -g "9-VARIANT_CALLING" ]] && [[ -g "9-VARIANT_CALLING/ANNOTATION" ]] ; then
	echo "Runing ANNOTATION"
	${S_PATH}/S9.5_vcf_results.sh
	${S_PATH}/S9.6_gtf_ref.sh
	Rscript ${S_PATH}/R9.9_Annotation
else
	if [[ -g "9-VARIANT_CALLING" ]] ; then
		echo "Creating ANNOTATION Directory"
		mkdir ANNOTATION
		echo "Runing ANNOTATION"
		${S_PATH}/S9.5_vcf_results.sh
		${S_PATH}/S9.6_gtf_ref.sh
		Rscript ${S_PATH}/R9.9_Annotation.R
	else
		echo "9-VARIANT_CALLING MISSING, copy directory from cluster after Variant Calling, variant calling script can be found in ${S_PATH}/S9_Variant_Calling.sh"
		exit -1
	fi

echo "STEP #9.5 ANNOTATION DONE"
echo "Final data table can be found in 9-VARIANT_CALLING/ANNOTATION/DATA_Variants_GATK.txt"


