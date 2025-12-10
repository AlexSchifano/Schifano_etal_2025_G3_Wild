setwd("~/Documents/C.elegans_sync/Final_files/01-DATA/9-VARIANT_CALLING/ANNOTATION/")

library(seqinr)
library(ade4)
library(adegenet)
library(tidyr)
library(stringr)
library(openxlsx)

Var1 = read.table("VCF_RES.txt",header = T)

## Removing mutations at pos 3229,3230, 3231, 3232 (this part was removed for ASR ==> always an insertion in VC)
remo_pos = c(3229,3230,3231,3232)
Var = Var1[Var1$POS_Full %in% remo_pos == F,]

## Add step to handle variants at same pos 
var_cor = which(grepl(",",Var$ALT))
for (f in var_cor) {
  print(f)
  print(Var[f,])
  other <- Var[f,]
  other1 <- gsub(",.*","",other)
  other2 <- gsub(".*,","",other)
  Var2 <- rbind(Var[1:f-1,],other1,other2,Var[seq(from = f+1, to = nrow(Var)),])
  Var <- Var2
}

Var$POS_Full <- as.numeric(Var$POS_Full)

## Only the CDS are protein coding, so I only keep the SNPs of CDS genes for annotation
co = colnames(Var)
CDS_Var = cbind(Var[Var$Gene_Type == "CDS",],0,0,0,0,0,0,0,0,0,0)
colnames(CDS_Var) = c(co,"POS_CDS","Effect_from","->","Effect_to","Synonymous","Transition","Radical","Triplet_Pos","4x_2x_Site","Codon")
CDS_Var$`->` <- "->"

## Amino acids groups (for labelling conservative vs radical)
{
  Aliphatic = c("G","A","V","L","I")
  Hydrocyl = c("S","C","U","T","M")
  Cyclic = c("P")
  Aromatic = c("F","Y","W")
  Basic = c("H","K","R")
  Acidic = c("D","E","N","Q")
  
  AA_groups = list(Aliphatic,Hydrocyl,Cyclic,Aromatic,Basic,Acidic)
}

## Annotation loop for CDS
for(R in 1:nrow(CDS_Var)) {
  
  ## Setting up the REF
  REF.file = paste("../REF/",CDS_Var[R,10],".fa",sep = "")
  REF.fa = read.fasta(REF.file,forceDNAtolower = F)
  REF=data.frame(seq(1:length(REF.fa[[1]])),REF.fa[[1]])
  colnames(REF) <- c("POS","SEQ")
  
  ## making the variant seq by replacing only the nuc in question
  VAR.fa = REF.fa
  VAR=data.frame(seq(1:length(VAR.fa[[1]])),VAR.fa[[1]])
  colnames(VAR) <- c("POS","SEQ")
  
  if (CDS_Var[R,8] == "SNP") {
    VAR[VAR$POS == CDS_Var[R,1],2] <- CDS_Var[R,4]
  }
  ## For INSERTIONS
  if (CDS_Var[R,8] == "INDEL" && nchar(CDS_Var[R,3]) < nchar(CDS_Var[R,4])) {
    
    ins_size = nchar(CDS_Var[R,4])-nchar(CDS_Var[R,3])
    
    INS <- data.frame(c(CDS_Var[R,1]:(CDS_Var[R,1]+ins_size)),str_split(CDS_Var[R,4],""))
    colnames(INS) <- c("POS","SEQ")
    
    VAR_tail <- data.frame(VAR[(CDS_Var[R,1]+ins_size):nrow(VAR),1]+ins_size, VAR[(CDS_Var[R,1]+ins_size):nrow(VAR),2])
    colnames(VAR_tail) <- c("POS","SEQ")
    
    VAR1 <- rbind(VAR[1:(CDS_Var[R,1]-1),],INS,VAR_tail)
    VAR <- VAR1
  }
  ## For DELETIONS ex:R=734
  if (CDS_Var[R,8] == "INDEL" && nchar(CDS_Var[R,3]) > nchar(CDS_Var[R,4])) {
    
    del_size = nchar(CDS_Var[R,3])-nchar(CDS_Var[R,4])
    
    VAR_tail <- data.frame(VAR[(CDS_Var[R,1]+del_size+1):nrow(VAR),1]-del_size, VAR[(CDS_Var[R,1]+del_size+1):nrow(VAR),2])
    colnames(VAR_tail) <- c("POS","SEQ")
    
    VAR1 <- rbind(VAR[1:(CDS_Var[R,1]),],VAR_tail)
    VAR <- VAR1
  }
  
  #### GTF STUFF
  {
  GTF.file = paste("../GTF/GTF_",CDS_Var[R,10],".fa.gtf",sep="")
  GTF = read.table(GTF.file)
  GTF[,6] = GTF[,3] - GTF[,2]+1
  colnames(GTF) = c("CHROM","START","END","GENE_ID","GENE_Type","SIZE")
  CDS_GTF = GTF[GTF$GENE_Type == "CDS",]
  
  P = seq(from = GTF$START[GTF$GENE_ID=="P"], to = GTF$END[GTF$GENE_ID=="P"])
  V = seq(from = GTF$START[GTF$GENE_ID=="V"], to = GTF$END[GTF$GENE_ID=="V"])
  ND6 = seq(from = GTF$START[GTF$GENE_ID=="ND6"], to = GTF$END[GTF$GENE_ID=="ND6"])
  ND4L = seq(from = GTF$START[GTF$GENE_ID=="ND4L"], to = GTF$END[GTF$GENE_ID=="ND4L"])
  W = seq(from = GTF$START[GTF$GENE_ID=="W"], to = GTF$END[GTF$GENE_ID=="W"])
  E = seq(from = GTF$START[GTF$GENE_ID=="E"], to = GTF$END[GTF$GENE_ID=="E"])
  rRNA12S = seq(from = GTF$START[GTF$GENE_ID=="12S_rRNA"], to = GTF$END[GTF$GENE_ID=="12S_rRNA"])
  S1 = seq(from = GTF$START[GTF$GENE_ID=="S1"], to = GTF$END[GTF$GENE_ID=="S1"])
  N = seq(from = GTF$START[GTF$GENE_ID=="N"], to = GTF$END[GTF$GENE_ID=="N"])
  Y = seq(from = GTF$START[GTF$GENE_ID=="Y"], to = GTF$END[GTF$GENE_ID=="Y"])
  ND1 = seq(from = GTF$START[GTF$GENE_ID=="ND1"], to = GTF$END[GTF$GENE_ID=="ND1"])
  ATP6 = seq(from = GTF$START[GTF$GENE_ID=="ATP6"], to = GTF$END[GTF$GENE_ID=="ATP6"])
  K = seq(from = GTF$START[GTF$GENE_ID=="K"], to = GTF$END[GTF$GENE_ID=="K"])
  L1 = seq(from = GTF$START[GTF$GENE_ID=="L1"], to = GTF$END[GTF$GENE_ID=="L1"])
  S2 = seq(from = GTF$START[GTF$GENE_ID=="S2"], to = GTF$END[GTF$GENE_ID=="S2"])
  ND2 = seq(from = GTF$START[GTF$GENE_ID=="ND2"], to = GTF$END[GTF$GENE_ID=="ND2"])
  I = seq(from = GTF$START[GTF$GENE_ID=="I"], to = GTF$END[GTF$GENE_ID=="I"])
  RR = seq(from = GTF$START[GTF$GENE_ID=="R"], to = GTF$END[GTF$GENE_ID=="R"])
  Q = seq(from = GTF$START[GTF$GENE_ID=="Q"], to = GTF$END[GTF$GENE_ID=="Q"])
  FF = seq(from = GTF$START[GTF$GENE_ID=="F"], to = GTF$END[GTF$GENE_ID=="F"])
  CTB1 = seq(from = GTF$START[GTF$GENE_ID=="CTB1"], to = GTF$END[GTF$GENE_ID=="CTB1"])
  L2 = seq(from = GTF$START[GTF$GENE_ID=="L2"], to = GTF$END[GTF$GENE_ID=="L2"])
  COX3 = seq(from = GTF$START[GTF$GENE_ID=="COX3"], to = GTF$END[GTF$GENE_ID=="COX3"])
  TT = seq(from = GTF$START[GTF$GENE_ID=="T"], to = GTF$END[GTF$GENE_ID=="T"])
  ND4 = seq(from = GTF$START[GTF$GENE_ID=="ND4"], to = GTF$END[GTF$GENE_ID=="ND4"])
  COX1 = seq(from = GTF$START[GTF$GENE_ID=="COX1"], to = GTF$END[GTF$GENE_ID=="COX1"])
  C = seq(from = GTF$START[GTF$GENE_ID=="C"], to = GTF$END[GTF$GENE_ID=="C"])
  M = seq(from = GTF$START[GTF$GENE_ID=="M"], to = GTF$END[GTF$GENE_ID=="M"])
  D = seq(from = GTF$START[GTF$GENE_ID=="D"], to = GTF$END[GTF$GENE_ID=="D"])
  G = seq(from = GTF$START[GTF$GENE_ID=="G"], to = GTF$END[GTF$GENE_ID=="G"])
  COX2 = seq(from = GTF$START[GTF$GENE_ID=="COX2"], to = GTF$END[GTF$GENE_ID=="COX2"])
  H = seq(from = GTF$START[GTF$GENE_ID=="H"], to = GTF$END[GTF$GENE_ID=="H"])
  rRNA16S = seq(from = GTF$START[GTF$GENE_ID=="16S_rRNA"], to = GTF$END[GTF$GENE_ID=="16S_rRNA"])
  ND3 = seq(from = GTF$START[GTF$GENE_ID=="ND3"], to = GTF$END[GTF$GENE_ID=="ND3"])
  ND5 = seq(from = GTF$START[GTF$GENE_ID=="ND5"], to = GTF$END[GTF$GENE_ID=="ND5"])
  A = seq(from = GTF$START[GTF$GENE_ID=="A"], to = GTF$END[GTF$GENE_ID=="A"])
  
  CDS.pos = c(ND6,ND4L,ND1,ATP6,ND2,CTB1,COX3,ND4,COX1,COX2,ND3,ND5)
  CDS.list = list(ND6,ND4L,ND1,ATP6,ND2,CTB1,COX3,ND4,COX1,COX2,ND3,ND5)
  names(CDS.list) <- c("ND6","ND4L","ND1","ATP6","ND2","CTB1","COX3","ND4","COX1","COX2","ND3","ND5")
  rm(ND6,ND4L,ND1,ATP6,ND2,CTB1,COX3,ND4,COX1,COX2,ND3,ND5)
  
  tRNA.pos = c(P,V,W,E,S1,N,Y,K,L1,S2,I,RR,Q,FF,L2,TT,C,M,D,G,H,A)
  tRNA.list = list(P,V,W,E,S1,N,Y,K,L1,S2,I,RR,Q,FF,L2,TT,C,M,D,G,H,A)
  names(tRNA.list) <- c("P","V","W","E","S1","N","Y","K","L1","S2","I","R","Q","FF","L2","TT","C","M","D","G","H","A")
  rm(P,V,W,E,S1,N,Y,K,L1,S2,I,RR,Q,FF,L2,TT,C,M,D,G,H,A)
  
  rRNA.pos = c(rRNA12S,rRNA16S)
  rRNA.list = list(rRNA12S,rRNA16S)
  names(rRNA.list) <- c("rRNA12s","rRNA16S")
  rm(rRNA12S,rRNA16S)
  
  Genes.list = c(CDS.list,tRNA.list,rRNA.list)
  }
  
  ## extract the CDS in REF
  REF_CDS = data.frame(REF[CDS.pos,])
  colnames(REF_CDS) <- c("POS","SEQ")
  
  ## extract the CDS in VAR
  VAR_CDS = data.frame(VAR[CDS.pos,])
  colnames(VAR_CDS) <- c("POS","SEQ")
  
  ## put CDS pos info in REF
  REF_CDS[,3] <- 1:nrow(REF_CDS)
  colnames(REF_CDS) <- c("POS","SEQ","CDS_POS")
  
  ## put CDS pos info in VAR
  VAR_CDS[,3] <- 1:nrow(VAR_CDS)
  colnames(VAR_CDS) <- c("POS","SEQ","CDS_POS")
  
  ## get the correction vector = vector of value to substract to the Full seq pos to get the pos from the CDS_string depending on the CDS
  FULL2CDS_POS = unique(REF_CDS[,1] - REF_CDS[,3])
  CDS_GTF = GTF[GTF$GENE_Type == "CDS",]
  CDS_GTF[,7] <- FULL2CDS_POS
  cog = colnames(GTF)
  colnames(CDS_GTF) = c(cog,"FULL_2_CDS_POS")
  
  ## put CDS pos correc
  gene = CDS_Var[R,2]
  corec = CDS_GTF[CDS_GTF$GENE_ID == gene,7]
  CDS_Var[R,11] <- CDS_Var[R,1] - corec 
  
  cds_pos = CDS_Var[R,11]
  cds_pos_AA =  ceiling(CDS_Var[R,11]/3)
  
  ## Adding triplet pos info ; 3rd pos /3 = full ; 2nd pos /3 = xx.66 ; 1st pos /3 = xx.33

    if (((cds_pos)/3) == ceiling(cds_pos/3)) { CDS_Var[R,18] <- 3 } 
    if (((cds_pos+1)/3) == ceiling(cds_pos/3)) { CDS_Var[R,18] <- 2 }
    if (((cds_pos+2)/3) == ceiling(cds_pos/3)) { CDS_Var[R,18] <- 1 }
  
  
  ## translate
  if (CDS_Var[R,8] == "SNP") {
  REF_AA = translate(REF_CDS$SEQ, numcode = 5)
  VAR_AA = translate(VAR_CDS$SEQ, numcode = 5)
  }
  
  ## Effect of CDS variant 
  ## Do something different for INDEL in CDS
  if (CDS_Var[R,8] == "SNP") {
    ref_aa <- REF_AA[cds_pos_AA]
    CDS_Var[R,12] <- ref_aa
    var_aa <- VAR_AA[cds_pos_AA] 
    CDS_Var[R,14] <- var_aa
  }
  
  ## For INDEL isolate gene to see if at the end of the beginning
  if (CDS_Var[R,8] == "INDEL") {
    gene_gtf <- GTF[GTF$GENE_ID == gene,]
    gene_pos = seq(from=gene_gtf$START,to=gene_gtf$END)
    
    REF_gene = REF[REF$POS %in% gene_pos,]
    
    ## For INS and DEL have to change the pos
    # DEL
    if (CDS_Var[R,8] == "INDEL" && nchar(CDS_Var[R,3]) > nchar(CDS_Var[R,4])) {
      gene_pos_del = seq(from=gene_gtf$START,to=(gene_gtf$END)-del_size)
      VAR_gene = VAR[VAR$POS %in% gene_pos_del,]
    }
    # INS
    if (CDS_Var[R,8] == "INDEL" && nchar(CDS_Var[R,3]) < nchar(CDS_Var[R,4])) {
      gene_pos_ins = seq(from=gene_gtf$START,to=(gene_gtf$END)+ins_size)
      VAR_gene = VAR[VAR$POS %in% gene_pos_ins,]
    }
    
    REF_gene_AA = translate(REF_gene$SEQ, numcode = 5)
    VAR_gene_AA = translate(VAR_gene$SEQ, numcode = 5)
    
    ## Correct to get pos in gene
    pos_gene <- CDS_Var[R,1] - gene_gtf$START + 1
    pos_gene_AA <- ceiling(pos_gene/3)
    
    ref_aas <- REF_gene_AA[(pos_gene_AA-1):(pos_gene_AA+15)]
    ref_aas <- ref_aas[which(ref_aas != "NA")]
    CDS_Var[R,12] <- paste(ref_aas,collapse = "")
    var_aas <- VAR_gene_AA[(pos_gene_AA-1):(pos_gene_AA+15)]
    var_aas <- var_aas[which(var_aas != "NA")]
    CDS_Var[R,14] <- paste(var_aas,collapse = "")

  }

  
  ## Synonymous or not
  if (CDS_Var[R,8] == "SNP") {
    if (ref_aa == var_aa) { CDS_Var[R,15] <- "Synonymous" } else { CDS_Var[R,15] <- "NON_Synonymous" }
  }
  if (CDS_Var[R,8] == "INDEL") {
    if (all(ref_aas == var_aas)) { CDS_Var[R,15] <- "Synonymous" } else { CDS_Var[R,15] <- "NON_Synonymous" }
  }
  
  ## Labeling Radical vs Conservative
  if (CDS_Var[R,8] == "SNP") {
    CDS_Var[R,17] <- "Radical"
    for (f in 1:6) {
      if (ref_aa %in% AA_groups[[f]] && var_aa %in% AA_groups[[f]]) { CDS_Var[R,17] <- "Conservative" }
    }
  }
  
  ## Labeling 4x or 2x site
  ## Here I only do it for the 3rd pos of codons
  ## 2x sites AAs are: F <-> L ; Y <-> * ; C <-> W ; H <-> Q ; I <-> M ; D <-> E ; N <-> K
  ## !!! Carefull, L is a special case, L has 6 codons, 4 start with C (4x site) and 2 start with T (2x sites)
  ## L is not in the next vector on purpose
    TwoXsites = c("F","Y","*","C","W","H","Q","I","M","D","E","N","K")
    
    if (CDS_Var[R,12] %in% TwoXsites && CDS_Var[R,18] == 3) { CDS_Var[R,19] <- "2x" }
    if (CDS_Var[R,12] == "L" && CDS_Var[R,18] == 3 && VAR_CDS$SEQ[cds_pos-2] == "T") { CDS_Var[R,19] <- "2x" }
    if (CDS_Var[R,12] == "L" && CDS_Var[R,18] == 3 && VAR_CDS$SEQ[cds_pos-2] == "C") { CDS_Var[R,19] <- "4x" }
    if ((CDS_Var[R,12] %in% TwoXsites) == F && CDS_Var[R,18] == 3 && CDS_Var[R,12] != "L") { CDS_Var[R,19] <- "4x" }
    if (CDS_Var[R,19] == 0) { CDS_Var[R,19] <- "N" }
  
  
  ## There is one case of 1st codon pos being 2x instead of 0x in the mt inv code
  ## Leucine (L) coded by TTA ; TTG ; CTA ; CTG   At those codons, first pos is 2x
  Codon_first_pos_exception = c("TTA","TTG","CTA","CTG")
  if (CDS_Var[R,18] == 1) {
    codon <- paste(VAR[CDS_Var[R,1],2],VAR[CDS_Var[R,1]+1,2],VAR[CDS_Var[R,1]+2,2],sep="")
    if (codon %in% Codon_first_pos_exception) { CDS_Var[R,19] <- "2x" } else { CDS_Var[R,19] <- "0x" }
  }
  
  ## At second codon pos, it is always 0x
  if (CDS_Var[R,18] == 2) { CDS_Var[R,19] <- "0x" }
  
  ## Add a column for codon in CDS
  if (CDS_Var[R,18] == 3) { CDS_Var[R,20] <- paste(VAR[CDS_Var[R,1]-2,2],VAR[CDS_Var[R,1]-1,2],VAR[CDS_Var[R,1],2],sep="") }
  if (CDS_Var[R,18] == 2) { CDS_Var[R,20] <- paste(VAR[CDS_Var[R,1]-1,2],VAR[CDS_Var[R,1],2],VAR[CDS_Var[R,1]+1,2],sep="") }
  if (CDS_Var[R,18] == 1) { CDS_Var[R,20] <- paste(VAR[CDS_Var[R,1],2],VAR[CDS_Var[R,1]+1,2],VAR[CDS_Var[R,1]+2,2],sep="") }

  
}

## Finishing the table
N_CDS_Var = cbind(Var[Var$Gene_Type != "CDS",],0,0,0,0,0,0,0,0,0,0)
colnames(N_CDS_Var) = c(co,"POS_CDS","Effect_from","->","Effect_to","Synonymous","Transition","Radical","Triplet_Pos","4x_2x_Site","Codon")

FINAL1 = rbind(N_CDS_Var,CDS_Var)

FINAL2 = unite(data = FINAL1, col = "MUTATION",c('REF','ALT'),sep = " -> ",remove = T)
FINAL3 = unite(data = FINAL2, col = "EFFECT",c('Effect_from','Effect_to'),sep = " -> ",remove = T)
FINAL = subset(FINAL3, select = -`->`)

## Transition vs transversion
Transi = c("A -> G","G -> A","C -> T","T -> C")

for (f in 1:nrow(FINAL)){
  if (FINAL[f,7] == "SNP"){
    if (FINAL[f,3] %in% Transi) { FINAL[f,13] <- "Transition" } else { FINAL[f,13] <- "Transversion" }
  }
}

write.xlsx(FINAL,file = "DATA_Variants_GATK.xlsx")
write.table(FINAL,file = "DATA_Variants_GATK.txt",sep="\t", quote = F, row.names = F)












