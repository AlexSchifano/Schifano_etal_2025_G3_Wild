setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")

###### Load the script with the functions I made ######
source("~/Documents/Scripts/R_FUNCTIONS.R")

library(openxlsx)
## Read the original file
Waneka = read.xlsx("Previous_results/Waneka/File_S1.xlsx")
library(stringr)

## Replace > by -> ; useful for mutation spectrum function
Waneka$Mutation <- str_replace(Waneka$`single.stranded.substitution.type.on.F-strand`, ">", " -> ")

## Making a better table
Waneka_New <- data.frame(Waneka$mtDNA.position,Waneka$Mutation,Waneka$reference,Waneka$variant,Waneka$genomic.region,Waneka$gene,Waneka$synonymous.or.nonysnonymous)
colnames(Waneka_New) <- c("POS","Mutation","REF","VAR","Gene_Type","Gene_ID","Effect")
Waneka_New$Effect[which(Waneka_New$Effect == "Syn ")] <- "Synonymous"
Waneka_New$Effect[which(Waneka_New$Effect == "Syn")] <- "Synonymous"

## Fix the Synonymous and Effect tagging
{
  library(tidyr)
  Eff_split <- separate(Waneka_New,Effect,into = c("Effect_from","Effect_to"),sep=">",remove = F)
  Eff_split$Synonymous <- 0
  Eff_split[is.na(Eff_split)] <- 0
  
  for (f in 1:nrow(Eff_split)) {
    if (Eff_split[f,"Effect_from"] == "Synonymous") { 
      Eff_split[f,"Synonymous"] <- "Synonymous"
      Eff_split[f,"Effect_to"] <- "Synonymous" 
    } else {
      if (Eff_split[f,"Effect_from"] != Eff_split[f,"Effect_to"]) { Eff_split[f,"Synonymous"] <- "NON_Synonymous" }
      if (Eff_split[f,"Effect_to"] == 0) { Eff_split[f,"Synonymous"] <- "." ; Eff_split[f,"Effect_to"] <- "." }
    }
  }
}

## Making the new table Eff_split but removing column 6
Waneka_New <- Eff_split[,-7]

## Add Transition info
{
  Ts.list <- c("C -> T","T -> C","A -> G","G -> A")
  Waneka_New$Transition <- 0
  for (f in 1:nrow(Waneka_New)){
    if (Waneka_New[f,"Mutation"] %in% Ts.list) { Waneka_New[f,"Transition"] <- "Transition" } else { Waneka_New[f,"Transition"] <- "Transversion" }
  }
}

## Radical vs Conservative
{
  ## Use to Eff_split_Syn table from earlier in th Waneka part of the script
  ## Amino acids groups (for labelling conservative vs radical)
  
  Aliphatic = c("G","A","V","L","I")
  Hydrocyl = c("S","C","U","T","M")
  Cyclic = c("P")
  Aromatic = c("F","Y","W")
  Basic = c("H","K","R")
  Acidic = c("D","E","N","Q")
  
  AA_groups = list(Aliphatic,Hydrocyl,Cyclic,Aromatic,Basic,Acidic)
  Waneka_New$Radical <- "."
  # I have to fix "Effect_to because there is a space
  Waneka_New$Effect_to <- str_remove(Waneka_New$Effect_to, pattern = " ")
  
  for (f in 1:nrow(Waneka_New)) {
    if (Waneka_New[f,"Synonymous"] == "NON_Synonymous") {
      Waneka_New[f,"Radical"] <- "Radical"
      for (g in 1:6) {
        if (Waneka_New[f,"Effect_from"] %in% AA_groups[[g]] && Waneka_New[f,"Effect_to"] %in% AA_groups[[g]]) { Waneka_New[f,"Radical"] <- "Conservative" }
        }
      }
    }
  }


## Get 2x and 4x info aswell as coding position
{
  # Split the codon change and effect
  Waneka_cod <- separate(Waneka,codon.change,into=c("Codon_from","Codon_to"),sep=">",remove=F)
  
  # Put position info and degeneracy
  Waneka_cod$Degeneracy <- "."
  Waneka_cod$Codon_position <- "."
  
  # Saving the colnames
  colna <- colnames(Waneka_New)
  
  # Binding the tables to have all the info in one
  Waneka_New <- cbind(Waneka_New, Waneka_cod$Codon_from, Waneka_cod$Codon_to, Waneka_cod$Degeneracy, Waneka_cod$Codon_position)
  
  # Fixing the colnames
  colnames(Waneka_New) <- c(colna,"Codon_from","Codon_to","Degeneracy","Codon_position")
  
  # Putting the actual info
  for (f in 1:nrow(Waneka_New)) {
    
    if (Waneka_New[f,"Gene_Type"]  == "CDS") {
      # get the codons (ref and mutant)
      cod_from <- str_split_1(Waneka_New[f,"Codon_from"],pattern = "")
      cod_to <- str_split_1(Waneka_New[f,"Codon_to"],pattern = "")
      
      # infer the codon pos of the mutation and store that data
      mut_cod_pos <- which(cod_from != cod_to)
      Waneka_New[f,"Codon_position"] <- mut_cod_pos
      
      # Add syn translation
      if (Waneka_New[f,"Effect_from"] == "Syn") { 
        Waneka_New[f,"Effect_from"] <- translate(cod_from,numcode = 5)
        Waneka_New[f,"Effect_to"] <- translate(cod_to,numcode = 5)
      }
      
      # Degeneracy
      degeneracy_tag_CDS(Sequence = cod_from, Codon_ = mut_cod_pos, Genetic_code = 5)
      Waneka_New[f,"Degeneracy"] <- Degeneracy
    } else {  }
    
  }
}

Waneka_New[11,7] <- "."
Waneka_New[is.na(Waneka_New)] <- "."
write.table(Waneka_New, file="Previous_results/Waneka_my_annotation.txt",sep="\t",quote = F,row.names = F,na = ".")


