setwd("~/Documents/C.elegans_sync/Final_files/01-DATA/8-COUNT_MUTATIONS/")
library(seqinr)
library(tidyr)


##### for testing #####

file="count_V301.fa"
anc = "V301"
comp = "V241"

#######################

args = commandArgs()
file = args[6]
anc = args[7]
comp = args[8]

cat("R inputs are",anc,comp,file)

## Reading the fasta with the 2 seq to compare
fasta = read.fasta(file = file, forceDNAtolower = F)

## Isolating Ancestor and Child (comp)
anc_fas1 = data.frame(seq(1:length(fasta[[anc]])),fasta[[anc]])
colnames(anc_fas1) = c("POS","SEQ")
comp_fas1 = data.frame(seq(1:length(fasta[[comp]])),fasta[[comp]]) 
colnames(comp_fas1) = c("POS","SEQ")

## The Alignment messed up and put a gap in the stop codon of ND6 instead of a few positions down
## Normal = TAA A ATT (ND6stop A ND4L start)
## Some lines have an insertion = TAA AA ATT
## But the alignment computes it as TAAAAATT
##                                  T-AAAATT
## I want to have                   TAA-AATT

## SO: If pos 548 = "-", I switch 548 with 551, now CDS ND6 is 115-549
## If pos 548 != "-", I do nothing and CDS is from 115-549

if (anc_fas1[548,2] == "-") {
  
  anc_fas <- anc_fas1
  anc_fas[548,2] <- anc_fas1[551,2]
  anc_fas[551,2] <- anc_fas1[548,2]
  
} else { anc_fas <- anc_fas1 }
if (comp_fas1[548,2] == "-") {
  
  comp_fas <- comp_fas1
  comp_fas[548,2] <- comp_fas1[551,2]
  comp_fas[551,2] <- comp_fas1[548,2]
  
} else { comp_fas <- comp_fas1 }


## Preparing the table of mutations
mut1=data.frame(matrix(ncol=5,nrow=1))
colnames(mut1) = c("POS","Mut_from","->","Mut_to","Mut_Type")

## Detecting the mutations
x=1
for (pos in 1:nrow(anc_fas)) {
  
  if (anc_fas[pos,2] != comp_fas[pos,2]) {
    mut1[x,1] <- pos
    mut1[x,2] <- anc_fas[pos,2]
    mut1[x,3] <- "->"
    mut1[x,4] <- comp_fas[pos,2]
    x=x+1
  }
  
}

if (is.na(mut1[1,1])) {
  q()
}
## Add mutation type info (SNP vs insertion vs deletion)
nuc = c("A","T","C","G")
for (m in 1:nrow(mut1)) {
  if (mut1[m,2] %in% nuc && mut1[m,4] %in% nuc) { mut1[m,5] <- "SNP" }
  if (mut1[m,2] == "-" && mut1[m,4] %in% nuc) { mut1[m,5] <- "INS" }
  if (mut1[m,2] %in% nuc && mut1[m,4] == "-") { mut1[m,5] <- "DEL" }
}

Mutation = unite(mut1[,2:4], "Mutation", remove = T, sep = " ")
mut = cbind(mut1[,1],Mutation,mut1[,5])
colnames(mut) = c("POS","Mutation","Mut_Type")

## GTF stuff
{
  
  GTF = read.table("GTF_blast.gtf",sep="\t")
  ## Gene size
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
  R = seq(from = GTF$START[GTF$GENE_ID=="R"], to = GTF$END[GTF$GENE_ID=="R"])
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
  
  tRNA.pos = c(P,V,W,E,S1,N,Y,K,L1,S2,I,R,Q,FF,L2,TT,C,M,D,G,H,A)
  tRNA.list = list(P,V,W,E,S1,N,Y,K,L1,S2,I,R,Q,FF,L2,TT,C,M,D,G,H,A)
  names(tRNA.list) <- c("P","V","W","E","S1","N","Y","K","L1","S2","I","R","Q","F","L2","T","C","M","D","G","H","A")
  rm(P,V,W,E,S1,N,Y,K,L1,S2,I,R,Q,FF,L2,TT,C,M,D,G,H,A)
  
  rRNA.pos = c(rRNA12S,rRNA16S)
  rRNA.list = list(rRNA12S,rRNA16S)
  names(rRNA.list) <- c("12S_rRNA","16S_rRNA")
  rm(rRNA12S,rRNA16S)
  
  Genes.list = c(CDS.list,tRNA.list,rRNA.list)

}

############### See effect of mut in CDS = translate

## extract the CDS
{
anc_CDS1 = data.frame(anc_fas[CDS.pos,])
colnames(anc_CDS1) = c("POS","SEQ")

comp_CDS1 = data.frame(comp_fas[CDS.pos,])
colnames(comp_CDS1) = c("POS","SEQ")
}
  
## remove the gaps, important for translation
{
anc_CDS = anc_CDS1[anc_CDS1$SEQ != "-",]
comp_CDS = comp_CDS1[comp_CDS1$SEQ != "-",]
}

## Correc pos after removing gaps by removing one to each pos after gap
{
  if (nrow(anc_CDS) != nrow(anc_CDS1)){
    gaps = anc_CDS1[anc_CDS1$SEQ == "-",]
    for (p in 1:nrow(gaps)) { 
      pos <- gaps[p,1]
      
      anc_CDS_pos <- c(anc_CDS[anc_CDS$POS < pos,1],anc_CDS[anc_CDS$POS > pos,1]-1)
      anc_CDS[,1] <- anc_CDS_pos
      }
  }
  
  if (nrow(comp_CDS) != nrow(comp_CDS1)){
    gaps = comp_CDS1[comp_CDS1$SEQ == "-",]
    for (p in 1:nrow(gaps)) { 
      pos <- gaps[p,1]
      
      comp_CDS_pos <- c(comp_CDS[comp_CDS$POS < pos,1],comp_CDS[comp_CDS$POS > pos,1]-1)
      comp_CDS[,1] <- comp_CDS_pos
    }
  }
}
  
## put CDS pos info
{
anc_CDS[,3] <- 1:nrow(anc_CDS)
colnames(anc_CDS) = c("POS","SEQ","CDS_POS")

comp_CDS[,3] <- 1:nrow(comp_CDS)
colnames(comp_CDS) = c("POS","SEQ","CDS_POS")
}

## get the correction vector = vector of value to substract to the Full seq pos to get the pos from the CDS_string depending on the CDS
{
FULL2CDS_POS = unique(anc_CDS[,1] - anc_CDS[,3])
CDS_GTF = GTF[GTF$GENE_Type == "CDS",]
CDS_GTF[,7] <- FULL2CDS_POS
co = colnames(GTF)
colnames(CDS_GTF) = c(co,"FULL_2_CDS_POS")
}

# translate
{
anc_AA = translate(anc_CDS$SEQ, numcode = 5)
comp_AA = translate(comp_CDS$SEQ,numcode = 5)
}

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

######### Make the final table
tab.mut = data.frame(matrix(nrow=0,ncol=16))
colnames(tab.mut) = c("POS_Full", "Gene_ID", "Mutation", "Gene_Type", "Variant_Freq", "Mut_Type", "Strain_Name", "Ancestor","Effect","Synonymous","Transition","Radical","Codon_Pos","Site_4x_2x","CDS_pos","Codon")

## ADD bunch on info to tab.mut
{
  transi = c("A -> G","G -> A","C -> T","T -> C")
  for (r in 1:nrow(mut)) {
    
    ## basic info we already had (Ancestor, Child, mutation, mutation type, position)
    {
    tab.mut[r,8] <- anc
    tab.mut[r,7] <- comp
    tab.mut[r,3] <- mut[r,2]
    tab.mut[r,1] <- mut[r,1]
    tab.mut[r,6] <- mut[r,3]
    }
    
    ## adding gene ID
    for (gene in 1:length(Genes.list)) {
      if (mut[r,1] %in% Genes.list[[gene]]) {
        tab.mut[r,2] <- names(Genes.list[gene])
      }
    }
    
    ## If it is not in any genes, it is in non coding region
    if (is.na(tab.mut[r,2])) { tab.mut[r,2] <- "NON_CODING" }
    
    if (tab.mut[r,2] %in% names(CDS.list)) { tab.mut[r,4] <- "CDS" }
    if (tab.mut[r,2] %in% names(tRNA.list)) { tab.mut[r,4] <- "tRNA" }
    if (tab.mut[r,2] %in% names(rRNA.list)) { tab.mut[r,4] <- "rRNA" }
    if (is.na(tab.mut[r,4])) { tab.mut[r,4] <- "NON_CODING" }
    
    ## Add transi ou transver to SNPs
    if (tab.mut[r,6] == "SNP" && tab.mut[r,3] %in% transi) { tab.mut[r,11] <- "Transition" }
    if (tab.mut[r,6] == "SNP" && (tab.mut[r,3] %in% transi) == F) { tab.mut[r,11] <- "Transversion" }
  }
}

for (r in 1:nrow(tab.mut)) {
  
  ## Put the CDS_POS info in tab.mut
  ## In english this loop is: 
  # we look at if a mutation is in a CDS, if it is, we look for the matching gene in CDS_GTF 
  # and substract the correction from the Full_Sea pos and store it in the CDS_POS (col 15 of tab.mut)
  # If the mutation doesn't occure in a CDS, no need for CDS pos since we won't check for translation effects
  if (tab.mut[r,4] == "CDS") {
    for (g in 1:nrow(CDS_GTF)) {
      if (tab.mut[r,2] == CDS_GTF[g,4]) { tab.mut[r,15] <- tab.mut[r,1] - CDS_GTF[g,7] }
    }
  }
  
  tab.mut[,15] <- as.numeric(tab.mut[,15])
  
  ## Using the same loop to add Codon pos
  if (tab.mut[r,4] == "CDS" && tab.mut[r,15]/3 == ceiling(tab.mut[r,15]/3)) { tab.mut[r,13] <- 3 }
  if (tab.mut[r,4] == "CDS" && (tab.mut[r,15]+1)/3 == ceiling(tab.mut[r,15]/3)) { tab.mut[r,13] <- 2 }
  if (tab.mut[r,4] == "CDS" && (tab.mut[r,15]+2)/3 == ceiling(tab.mut[r,15]/3)) { tab.mut[r,13] <- 1 }
  
  ## Same loop for evrything that has to do with AA (Effect, Syn, Radical, 4x site...)
  if (tab.mut[r,4] == "CDS") {
    
    # fetching the AA at mut pos
    pos_AA = ceiling(as.numeric(tab.mut[r,15])/3)
    Eff_anc = anc_AA[[pos_AA]]
    Eff_comp = comp_AA[[pos_AA]]
    
    # putting Effect in tab.mut
    tab.mut[r,9] <- paste(Eff_anc, "->", Eff_comp)
    
    # putting if Synonymous or not
    if (Eff_anc == Eff_comp) { tab.mut[r,10] <- "Synonymous" } else { tab.mut[r,10] <- "NON_Synonymous" }
    
    # radical or conservative change
    # if Synonymous = automatically conservative
    # if not check if in the same AA groups, if not => Radical
    if (tab.mut[r,10] == "Synonymous") { tab.mut[r,12] <- "Conservative" } else {
      for (g in 1:length(AA_groups)) {
        if (Eff_anc %in% AA_groups[[g]] && Eff_comp %in% AA_groups[[g]]) { tab.mut[r,12] <- "Conservative" }
      }
      if (is.na(tab.mut[r,12])) { tab.mut[r,12] <- "Radical" }
    }
    
    ## Labeling 4x or 2x site
    ## Here I only do it for the 3rd pos of codons
    ## 2x sites AAs are: F <-> L ; Y <-> * ; C <-> W ; H <-> Q ; I <-> M ; D <-> E ; N <-> K
    ## !!! Carefull, L is a special case, L has 6 codons, 4 start with C (4x site) and 2 start with T (2x sites)
    ## L is not in the next vector on purpose
    TwoXsites = c("F","Y","*","C","W","H","Q","I","M","D","E","N","K")
    
    if (Eff_anc %in% TwoXsites && tab.mut[r,13] == 3) { tab.mut[r,14] <- "2x" }
    if (Eff_anc == "L" && tab.mut[r,13] == 3 && anc_fas[tab.mut[r,1]-2,2] == "T") { tab.mut[r,14] <- "2x" }
    if (Eff_anc == "L" && tab.mut[r,13] == 3 && anc_fas[tab.mut[r,1]-2,2] == "C") { tab.mut[r,14] <- "4x" }
    if ((Eff_anc %in% TwoXsites) == F && tab.mut[r,13] == 3 && Eff_anc != "L") { tab.mut[r,14] <- "4x" }
    
    ## There is one case of 1st codon pos being 2x instead of 0x in the mt inv code
    ## Leucine (L) coded by TTA ; TTG ; CTA ; CTG   At those codons, first pos is 2x
    Codon_first_pos_exception = c("TTA","TTG","CTA","CTG")
    if (tab.mut[r,13] == 1) {
      codon <- paste(comp_fas[tab.mut[r,1],2],comp_fas[tab.mut[r,1]+1,2],comp_fas[tab.mut[r,1]+2,2],sep="")
      if (codon %in% Codon_first_pos_exception) { tab.mut[r,14] <- "2x" } else { tab.mut[r,14] <- "0x" }
    }
    
    ## At second codon pos, it is always 0x
    if (tab.mut[r,13] == 2) { tab.mut[r,14] <- "0x" }
   
  }
  
  ## Add the codon for CDS mutation
  if (tab.mut[r,4] == "CDS") {
    if (tab.mut[r,13] == 1) { tab.mut[r,"Codon"] <- paste(comp_fas[tab.mut[r,1],2],comp_fas[tab.mut[r,1]+1,2],comp_fas[tab.mut[r,1]+2,2],sep="") }
    if (tab.mut[r,13] == 2) { tab.mut[r,"Codon"] <- paste(comp_fas[tab.mut[r,1]-1,2],comp_fas[tab.mut[r,1],2],comp_fas[tab.mut[r,1]+1,2],sep="") }
    if (tab.mut[r,13] == 3) { tab.mut[r,"Codon"] <- paste(comp_fas[tab.mut[r,1]-2,2],comp_fas[tab.mut[r,1]-1,2],comp_fas[tab.mut[r,1],2],sep="") }
  }
  
  
}

write.table(tab.mut, file = "Count_mut_results.txt", sep = "\t", quote = F, row.names = F, col.names = F, append = T)
#write.csv(tab.mut, file = "Mutations_Tree_ASR.xlsx", quote = F, row.names = F, col.names = F, append = T)
