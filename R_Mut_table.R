setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")

###### Reading the various mutation lists and merging it ###### 
{
  Var = read.table("../01-DATA/9-VARIANT_CALLING/ANNOTATION/DATA_Variants_GATK_REF_POS.txt",header = T, sep="\t")
  ASR = read.table("../01-DATA/8-COUNT_MUTATIONS/Count_Mut_results_REF_POS.txt", header = T, sep="\t")
  Table_VC = read.table("../01-DATA/9-VARIANT_CALLING/Table_VC_All.txt")
  MutServe = read.table("../01-DATA/mutserve/ANNOTATION/DATA_Variants_mutserve.txt",header = T, sep="\t")
  
  ## Removing some un-usefull columns
  Var = subset(Var, select = -`QUAL`)
  Var = subset(Var, select = -`POS_CDS`)
  ASR = subset(ASR, select = -`CDS_pos`)
  
  ## Synchronizing colnames
  co = colnames(ASR)
  colnames(Var) <- co
  
  ## Making a MEGA table
  MEGA = rbind(Var,ASR)
  MEGA[is.na(MEGA)] <- "N"
}

##### Fixing indel positions #####
{
  MEGA_INDEL <- MEGA[MEGA$Mut_Type != "SNP",]
  
  # I only need to fix the position of the indels called within the tree, they are called INS or DEL (not INDEL)
  MEGA_INDEL_fix <- MEGA_INDEL[MEGA_INDEL$Mut_Type != "INDEL",]
  table(MEGA_INDEL_fix$POS_Full)
  ### There are few enough positions for me to fix it by hand...
  ## Pos 18 is a homopolymer, we will make it fit the GATK format
  {
    D18 <- MEGA_INDEL_fix[MEGA_INDEL_fix$POS_Full == 18,]
    D18 <- D18[D18$Mutation == "A -> -",]
    D18$POS_Full <- 16
    D18$POS_REF <- 16
    D18$Mutation <- "TA -> T"
    
    I18 <- MEGA_INDEL_fix[MEGA_INDEL_fix$POS_Full == 18,]
    I18 <- I18[I18$Mutation == "- -> A",]
    I18$POS_Full <- 16
    I18$POS_REF <- 16
    I18$Mutation <- "T -> TA"
  }
  ## Pos 551
  {
    D551 <- MEGA_INDEL_fix[MEGA_INDEL_fix$POS_Full == 551,]
    D551 <- D551[D551$Mutation == "A -> -",]
    D551$POS_Full <- 549
    D551$POS_REF <- 546
    D551$Mutation <- "TA -> T"
    
    I551 <- MEGA_INDEL_fix[MEGA_INDEL_fix$POS_Full == 551,]
    I551 <- I551[I551$Mutation == "- -> A",]
    I551$POS_Full <- 549
    I551$POS_REF <- 546
    I551$Mutation <- "T -> TA"
  }
  ## Pos 804
  {
    D804 <- MEGA_INDEL_fix[MEGA_INDEL_fix$POS_Full == 804,]
    D804$POS_Full <- 803
    D804$POS_REF <- 800
    D804$Mutation <- "GT -> G"
  }
  ## Pos 894
  {
    D894 <- MEGA_INDEL_fix[MEGA_INDEL_fix$POS_Full == 894,]
    D894$POS_Full <- 888
    D894$POS_REF <- 885
    D894$Mutation <- "GT -> T"
  }
  ## Pos 1755
  {
    D1755 <- MEGA_INDEL_fix[MEGA_INDEL_fix$POS_Full == 1755,]
    D1755 <- D1755[D1755$Mutation == "T -> -",]
    D1755$POS_Full <- 1754
    D1755$POS_REF <- 1750
    D1755$Mutation <- "AT -> A"
    
    I1755 <- MEGA_INDEL_fix[MEGA_INDEL_fix$POS_Full == 1755,]
    I1755 <- I1755[I1755$Mutation == "- -> T",]
    I1755$POS_Full <- 1754
    I1755$POS_REF <- 1750
    I1755$Mutation <- "A -> AT"
  }
  ## Pos 1756
  {
    D1756 <- MEGA_INDEL_fix[MEGA_INDEL_fix$POS_Full == 1756,]
    D1756$POS_Full <- 1754
    D1756$POS_REF <- 1750
    D1756$Mutation <- "AT -> A"
  }
  ## Pos 7730
  {
    D7730 <- MEGA_INDEL_fix[MEGA_INDEL_fix$POS_Full == 7730,]
    D7730$POS_Full <- 7729
    D7730$POS_REF <- 7736
    D7730$Mutation <- "TA -> T"
  }
  ## Pos 9577
  {
    D9577 <- MEGA_INDEL_fix[MEGA_INDEL_fix$POS_Full == 9577,]
    D9577$POS_Full <- 9576
    D9577$POS_REF <- 9582
    D9577$Mutation <- "T -> TA"
  }
  ## Pos 9632
  {
    D9632 <- MEGA_INDEL_fix[MEGA_INDEL_fix$POS_Full == 9632,]
    D9632$POS_Full <- 9631
    D9632$POS_REF <- 9636
    D9632$Mutation <- "G -> GT"
  }
  
  ## Re merge the indel table
  MEGA_INDEL_nofix <- MEGA_INDEL[MEGA_INDEL$Mut_Type == "INDEL",]
  MEGA_INDEL <- rbind(MEGA_INDEL_nofix,D18,I18,D551,I551,D804,D894,D1755,I1755,D1756,D7730,D9577,D9632)
  
  ## Re merge the whole table
  MEGA_snp <- MEGA[MEGA$Mut_Type == "SNP",]
  MEGA <- rbind(MEGA_snp,MEGA_INDEL)
}

###### Removing parallel mutations that occurred in same mitotype lines ###### 
{
  library(openxlsx)
  ## Here I am listing the parallel mutations to check if we still have contamination even after removing the contaminated lines
  ## I also look at the mutations that are shared within a mitotype in order to only account for one set of those
  ASR = read.table("../01-DATA/8-COUNT_MUTATIONS/Count_mut_results.txt", header = T, sep="\t")
  Shared <- data.frame(table(MEGA$POS_Full))
  write.table(Shared,file = "shared.txt",quote = F,row.names = F)
  # Shared.txt has the mutation count for each position
  Shared <- read.table("shared.txt",header = T)
  file.remove("shared.txt")
  # subset it to only take positions found more than once
  Shared_2p <- Shared[Shared$Freq > 1,]
  ## Loop to get parallel mutations
  Paralel_muts <- data.frame()
  for (pos in Shared_2p$Var1) {
    mega <- MEGA[MEGA$POS_Full == pos,]
    for(i in 1:nrow(mega)) {
      for (f in 1:nrow(mega)) {
        if (i != f) { if(sum(mega[i,c(1:3,8)] == mega[f,c(1:3,8)]) == 4) {
          Paralel_muts <- rbind(Paralel_muts,mega[f,])
        }}}}}
  # removing natural repeats coming from the loop
  Paralel_muts <- unique(Paralel_muts)
  # Count for each position in Paralel_muts
  Strain_shared_mut <- data.frame(table(Paralel_muts$Strain_Name))
  colnames(Strain_shared_mut) <- c("Strain","NB_Shared_Mut")
  write.table(Strain_shared_mut,file="List_strains_with_shared_mutations.txt",sep="\t",quote = F,row.names = F)
  Strain_shared_mut <- read.table("List_strains_with_shared_mutations.txt",sep="\t",header = T)
  file.remove("List_strains_with_shared_mutations.txt")
  # Subset Paralel_muts to get the exact number of paralele sites
  Paralel_muts_sub <- Paralel_muts[,1:3]
  Paralel_muts_sub <- unique(Paralel_muts_sub)
  # ==> There are 522 positions involved
  ## Loop to add variant freq info
  MEGA_shared <- MEGA[MEGA$POS_Full %in% Paralel_muts$POS_Full,]
  Strain_shared_mut$Min_VAF <- 0
  Strain_shared_mut$Mean_VAF <- 0
  Strain_shared_mut$Max_VAF <- 0
  row.names(Strain_shared_mut) <- Strain_shared_mut$Strain
  for (strain in Strain_shared_mut$Strain) {
    Strain_shared_mut[strain,"Mean_VAF"] <- mean(as.numeric(MEGA_shared$Variant_Freq[MEGA_shared$Strain_Name == strain]))
    Strain_shared_mut[strain,"Min_VAF"] <- min(as.numeric(MEGA_shared$Variant_Freq[MEGA_shared$Strain_Name == strain]))
    Strain_shared_mut[strain,"Max_VAF"] <- max(as.numeric(MEGA_shared$Variant_Freq[MEGA_shared$Strain_Name == strain]))
  }

  ## Only keeping the strains and no nodes
  Strain_shared_mut_str <- Strain_shared_mut[Strain_shared_mut$Strain %in% ASR$Strain_Name == F,]
    
  ## Only keeping strains with 2 or more shared mutations
  Strain_shared_mut_2plus <- Strain_shared_mut[Strain_shared_mut$NB_Shared_Mut > 1,]
  
  # Adding ancestor info since some parallel mutations come from strains of the same cluster
  Strain_shared_mut_str$Ancestor <- ""
  Table_VC <- read.table("../01-DATA/7-ARPIP/Output/Table_VC_ALL.txt")
  for (f in 1:nrow(Strain_shared_mut_str)) {
    str <- Strain_shared_mut_str[f,"Strain"]
    Strain_shared_mut_str[f,"Ancestor"] <- Table_VC$V2[Table_VC$V1 == str]
  }
  # Counts of shared ancestor to remove parallels that escaped the method
  Anc_count_shared <- data.frame(table(Strain_shared_mut_str$Ancestor))
  # Remove the ones with only 1 occurrence (no fake parallelism)
  Anc_paralel <- Anc_count_shared[Anc_count_shared$Freq > 1,]
  # Get a table with the concerned strains
  Strains_fakePara <- Strain_shared_mut_str[Strain_shared_mut_str$Ancestor %in% Anc_paralel$Var1,]
  # I now make a list of the lines I want to keep based on the ones that cluster together in the same mitotype (and have parallel mutations)
  # I keep one line per mitotype
  # I do it by hand in order to not remove lines with potential unique variants at lower frequencies
  # I looked at the table Strains_fakePara and sorted by Ancestor. When several strains shared the same ancestor and had similar counts  
  # and VAF I picked one of them to keep, I confirmed those by looking at the Paralel_muts table and cheking if the mutations were trully parallel
  # in some cases the mutations were not parallel within a mitotype but with other lines without any pattern. 
  list_strains <- c("ECA2568","NIC207","ED3005","CX11315","JU3144","ECA348","JU792","ECA640","QX1792",
                    "ECA1269","ECA760","ECA1977","NIC1806","JU4082","NIC256","CB4932","JU1934","JU1896",
                    "NIC261","JU1246","NIC268","ECA1887","ECA1257","ECA1693","ECA723","ECA1237","ECA1252",
                    "ECA1232","ECA363","QX1794","ECA1805","ECA2367","ECA2043","QX1211","ECA1717","ECA1851","ED3049")
  ## Now I make the list of strains to remove the mutations of
  list_remove <- Strains_fakePara$Strain[(Strains_fakePara$Strain %in% list_strains) == F]
  ## Now the list of mutations that are parallel in those lines so that I don't remove the non-parallel ones
  Paralel_muts_rm <- Paralel_muts[Paralel_muts$Strain_Name %in% list_remove,]
  ## Removing those mutations from MEGA
  library(dplyr)
  MEGA_fixed <- anti_join(MEGA,Paralel_muts_rm)
  if(nrow(MEGA) - nrow(Paralel_muts_rm) == nrow(MEGA_fixed)) {
    write.table(MEGA_fixed,sep="\t",col.names = T,row.names = F,file = "R_txt_files/Mutation_List.txt",quote = F)
    print("Mutation list without within mitotype parallelisms writen in Mutation_List.txt in folder R_txt_files")
    MEGA <- MEGA_fixed
    }else{print("ERROR in the process, more mutations removed than should be")}
}

###### Making sub tables #####
{
  # I rename the first column
  cols <- colnames(MEGA)
  cols[1] <- "POS_Anc"
  colnames(MEGA) <- cols
  
  MEGA_CDS = MEGA[MEGA$Gene_Type == "CDS",]
  MEGA_tRNA = MEGA[MEGA$Gene_Type == "tRNA",]
  MEGA_rRNA = MEGA[MEGA$Gene_Type == "rRNA",]
  MEGA_NON_CODING = MEGA[MEGA$Gene_Type == "NON_CODING",]
  MEGA_Fourfold = MEGA[MEGA$Site_4x_2x == "4x",]
  MEGA_Twofold = MEGA[MEGA$Site_4x_2x == "2x",]
  MEGA_NonDegenerate = MEGA[MEGA$Site_4x_2x == "0x",]
  MEGA_INDEL = MEGA[MEGA$Mut_Type != "SNP",]
  MEGA_Var <- MEGA[MEGA$Strain_Name %in% Var$Strain_Name,]
}

###### Adding a column for GENE_POS (= FULL_POS - GENE_START + 1) to MEGA_CDS ###### 
{
  V170_REF_GTF = read.table("../01-DATA/9-VARIANT_CALLING/GTF/GTF_V170.fa.gtf")
  V170_REF_GTF <- subset(V170_REF_GTF, select = -V1)
  V170_REF_GTF[,5] <- V170_REF_GTF[,2] - V170_REF_GTF[,1] +1
  colnames(V170_REF_GTF) = c("START","END","Gene_ID","Gene_Type","SIZE")
  V170_CDS <- V170_REF_GTF[V170_REF_GTF$Gene_Type == "CDS",]
  
  MEGA_CDS$POS_GENE <- MEGA_CDS$POS_Full
  for (gene in V170_CDS$Gene_ID) {
    MEGA_CDS$POS_GENE[MEGA_CDS$Gene_ID == gene] <- MEGA_CDS$POS_Full[MEGA_CDS$Gene_ID == gene] - V170_CDS$START[V170_CDS$Gene_ID == gene] + 1}
}

##### Writing the subseted tables #####
{
  write.table(MEGA,sep="\t",col.names = T,row.names = F,file = "R_txt_files/Mutation_List.txt",quote = F)
  write.xlsx(MEGA, file = "R_txt_files/Mutation_List.xlsx")
  write.table(MEGA_CDS,sep="\t",col.names = T,row.names = F,file = "R_txt_files/Mutation_CDS.txt",quote = F)
  write.table(MEGA_tRNA,sep="\t",col.names = T,row.names = F,file = "R_txt_files/Mutation_tRNA.txt",quote = F)
  write.table(MEGA_rRNA,sep="\t",col.names = T,row.names = F,file = "R_txt_files/Mutation_rRNA.txt",quote = F)
  write.table(MEGA_NON_CODING,sep="\t",col.names = T,row.names = F,file = "R_txt_files/Mutation_NON_CODING.txt",quote = F)
  write.table(MEGA_Fourfold,sep="\t",col.names = T,row.names = F,file = "R_txt_files/Mutation_Fourfold.txt",quote = F)
  write.table(MEGA_Twofold,sep="\t",col.names = T,row.names = F,file = "R_txt_files/Mutation_Twofold.txt",quote = F)
  write.table(MEGA_NonDegenerate,sep="\t",col.names = T,row.names = F,file = "R_txt_files/Mutation_NonDegenerate.txt",quote = F)
  write.table(MEGA_Var, sep="\t",col.names = T,row.names = F,file = "R_txt_files/Variants_List.txt",quote = F)
}
