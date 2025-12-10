setwd("~/Documents/C.elegans_sync/Final_files/01-DATA/9-VARIANT_CALLING/ANNOTATION/")

library(vcfR)
library(ape)
library(tidyr)
library(stringr)

## Script to read vcf files and only keep the SNPs and put everything in a table with the strain data
## Script used with bash script Ssnp_vcf_strains.sh
## Further data processing in script R_Anotation.R
## Script looped with VCF names

## Testing ##########
vcf.name = "../VCF/AB1.vcf"
Strain.name = "AB1"
gtf = "../GTF/GTF_V307.fa.gtf"
anc = "V2"
#####################

args = commandArgs()
vcf.name = args[6]
Strain.name = args[7]
gtf = args[8]
anc = args[9]

print(Strain.name)

gtf=read.table(gtf)
colnames(gtf) <- c("CHROM","START","END","ID","TYPE")
vcf = read.vcfR(file = vcf.name)


###  Create a chromR object  ###
# name : Chrom name for plotting  # seq : sequence as DNAbin object (see adegenet package)
mtdna <- create.chromR(name = "plop", vcf=vcf)
print(Strain.name)
# spliting the FORMAT column
muts.0 = data.frame(mtdna@vcf@fix,0,vcf@gt[,2])
colnames(muts.0) = c(colnames(mtdna@vcf@fix),"Gene_Type","FORMAT")
FOR = c("GT","AD","FORMAT_DP","GQ","PL","VAF")
INF = c("AC","AF","AN","BaseRankSum","DP","ExcessHet","FS","MLEAC","MLEAF","MQ","MQRankSum","QD","ReadPosRankSum","SOR","TYPE")

muts.1 = separate(data = muts.0, col = FORMAT, into = FOR, sep = ":")

## !!!! CAREFULL Only run this loop once !!!! ##
## Loop to reorder the columns because different number of arguments in the INFO col depending on variant

## Run this to reset muts

I = as.data.frame(muts.1$INFO)
colnames(I) = "INFO"
SEP = separate(data = I, col = INFO,into = INF, sep = ";")
TT = c("TYPE=INDEL","TYPE=SNP")
for (r in seq(1:nrow(SEP))) {
  for (c in seq(1:ncol(SEP))) {
    if (SEP[r,c] %in% TT) { 
      SEP[r,c]
      SEP[r,15] <- SEP[r,c] }
  }
}

## Only keep the TYPE from INFO
muts <- muts.1
muts$INFO <- SEP$TYPE

## Removing the XX= in each line

muts[,8] <- str_remove(muts[,8], "TYPE=")

muts[is.na(muts)] <- "NN"

### Loop to ad gene ID and Gene_Type to the mut table

muts$POS <- as.numeric(muts$POS)
CDS.lst = gtf[gtf$TYPE == "CDS",4]
tRNA.lst = gtf[gtf$TYPE == "tRNA",4]
rRNA.lst = gtf[gtf$TYPE == "rRNA",4]
GENE.lst = c(CDS.lst,tRNA.lst,rRNA.lst)

for (f in seq(1:nrow(muts))) {
  
  for (x in seq(1:nrow(gtf))) {
    
    START = as.numeric(gtf[x,2])
    END = as.numeric(gtf[x,3])
    
    if (muts[f,2] %in% seq(from = START, to = END)) {
      
      muts[f,3] <- gtf[x,4]
    }
  }
  
  if (muts[f,3] == "NN") { muts[f,3] <- "NON_CODING" }
}

for (g in seq(1:nrow(muts))) {
  
  if (muts[g,3] %in% CDS.lst) {
    muts[g,9] <- "CDS" }
  
  if (muts[g,3] %in% tRNA.lst) {
    muts[g,9] <- "tRNA" }
  
  if (muts[g,3] %in% rRNA.lst) {
    muts[g,9] <- "rRNA" }
  
  if (muts[g,3] %in% GENE.lst == F) { 
    muts[g,9] <- "NON_CODING" }
}

## Add Strain and Anc info
muts[,16] <- Strain.name
muts[,17] <- anc

## restructuring the table
# "POS_Full Gene_ID REF ALT QUAL Gene_Type Variant_Freq Mut_Type Strain_Name Ancestor"
muts.fin = cbind(muts[,2:6],muts[,9],muts[,15],muts[,8],muts[,16:17])

## Write table
write.table(muts.fin, file="VCF_RES.txt",append=TRUE,quote=F,row.names = F,sep = "\t",col.names = F)

















