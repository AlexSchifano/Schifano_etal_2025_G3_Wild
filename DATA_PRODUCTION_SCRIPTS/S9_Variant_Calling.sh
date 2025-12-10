#!/bin/bash -l

#SBATCH -A naiss2023-5-394

#SBATCH -p core

#SBATCH -n 20

#SBATCH -t 10:00:00

#SBATCH -J iqtree_elegans_FULL

set -euo pipefail 

## Variant Calling pipeline for calling variants comparing to ASR sequences
## Before running make sure there are folders called "REF" "BAM" "VCF" ==> mkdir REF BAM VCF

module load bioinfo-tools
module load FastQC/0.11.9
module load bwa/0.7.17
module load samtools/1.17
module load GATK/4.3.0.0
module load bcftools/1.17
module load python/3.9.5
module load SeqKit/2.4.0

cd /proj/snic2022-6-224/Alex//Final_files/01-DATA/9-VARIANT_CALLING

# copy all the seqs
cp ../7-ARPIP/Output/Anc.fasta .

# make a version without gaps to use as reference
sed 's/-//g' Anc.fasta > No_GAP_Anc.fasta

while IFS=' ' read comp anc
do       
  NODE=${anc}
  LEAF=${comp}
  REF=${anc}.fa
  BAM=${comp}_MtDNA.bam



########### REFERENCE ############

# get the ancestral node seq that will be used as a reference for the variant calling
# Need to make a file with the name of seq to grep from msa
echo "$anc" > extract.txt
seqkit grep -n -f extract.txt No_GAP_Anc.fasta -o $REF
rm extract.txt

# remove the gaps in the ref and change the name
sed '1 s/^.*$/>MtDNA/' $REF > REF/${REF}
rm $REF

# Index the reference genome (NodeX.fa)
bwa index -a bwtsw REF/$REF
samtools faidx REF/$REF

if [[ -f REF/${anc}.dict ]];
then
  echo "Exists"
else
  gatk --java-options -Xmx7g CreateSequenceDictionary -R REF/$REF -O REF/${anc}.dict
fi

############### REFERENCE READY ##################


############### BAM ##############

# get the bam files that will be used for the variant calling, one for each descendent of the node, comp1 et comp2
cp ../1-MtDNA_BAM/$BAM .

# index the LEAF bam
samtools index $BAM

# remove Dloop
samtools view -b $BAM "MtDNA:1-13328" > ${LEAF}_trim1.bam
samtools index ${LEAF}_trim1.bam


## Need to re-align the reads to the new REF

# convert the bam back to a fastq
samtools fastq ${LEAF}_trim1.bam > ${comp}.fq

# remove the original BAM and all the indexes
rm $BAM
rm *.bam*

# Align the reads to the reference
bwa mem -R "@RG\\tID:ALN2_${comp}\\tSM:${comp}\\tLB:illumina\\tPL:illumina" -t 1 REF/$REF ${LEAF}.fq | samtools sort > BAM/$BAM

# Index the new BAM
samtools index BAM/$BAM

############## BAM FILE READY ################

# MarkDuplicates
gatk --java-options -Xmx7g MarkDuplicates -I BAM/$BAM -O BAM/${LEAF}_marked.bam -M DuplicateMetrics/${LEAF}_marked_dup_metrics.txt

# Index the marked BAM
samtools index BAM/${LEAF}_marked.bam

##################################################


############### VARIANT CALLING AND ADD INFO ##############

# Variant calling with adding info about VAF and Mut Type
gatk --java-options -Xmx7g HaplotypeCaller -R REF/$REF -I BAM/${LEAF}_marked.bam -O temp_${LEAF}.vcf
bcftools +fill-tags temp_${LEAF}.vcf -o temp_${LEAF}_2.vcf -- -t FORMAT/VAF
bcftools +fill-tags temp_${LEAF}_2.vcf -o VCF/${LEAF}.vcf -- -t INFO/TYPE

rm ${LEAF}*
rm temp_*

done < Table_VC_All.txt

echo "All vcf can be found in /proj/snic2022-6-224/Alex/C.elegans_sync/Mutations/FullSeq_IQtree_ASR/WASR/VariantCalling/VCF"
