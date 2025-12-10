#!/bin/bash

## Run S9.6 first (I know it is not in order sorry)
## Script to extract mut info from vcf ouptu from GATK Haplotype Caller

set -euo pipefail

cd /Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/9-VARIANT_CALLING/ANNOTATION

ls ../VCF > VCF.list

if [[ -f VCF_RES ]] ; then
	rm VCF_RES.txt
fi

echo "POS_Full Gene_ID REF ALT QUAL Gene_Type Variant_Freq Mut_Type Strain_Name Ancestor" | sed "s/ /\t/g" > VCF_RES.txt

for f in $(cat VCF.list) ; do
	
	strain=$(basename $f .vcf)
	vcf=../VCF/${f}
	anc=$(grep -w "$strain" ../Table_VC_All.txt | sed "s/$strain\t//g")
	gtf=../GTF/GTF_${anc}.fa.gtf

	Rscript ../../0-Scripts/R9.5_vcf_results.R $vcf $strain $gtf $anc ; 
done
