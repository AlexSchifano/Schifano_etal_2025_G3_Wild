setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")

###### Load the script with the functions I made ######
source("~/Documents/Scripts/R_FUNCTIONS.R")

###### Load the mutation lists ######
{
  MEGA <- read.table("R_txt_files/Mutation_List.txt",header = T,sep="\t")
  MEGA_CDS <- read.table("R_txt_files/Mutation_CDS.txt",header = T,sep="\t")
  MEGA_tRNA <- read.table("R_txt_files/Mutation_tRNA.txt",header = T,sep="\t")
  MEGA_rRNA <- read.table("R_txt_files/Mutation_rRNA.txt",header = T,sep="\t")
  MEGA_NON_CODING <- read.table("R_txt_files/Mutation_NON_CODING.txt",header = T,sep="\t")
  MEGA_0x <- read.table("R_txt_files/Mutation_NonDegenerate.txt",sep="\t",header=T)
  MEGA_2x <- read.table("R_txt_files/Mutation_Twofold.txt",sep="\t",header=T)
  MEGA_4x <- read.table("R_txt_files/Mutation_Fourfold.txt",sep="\t",header=T)
}

##### Mutation stats #####
{
  ##### Mutation proportions #####
  {
    ##### Per Codon position #####
    {
      ## Counting number of mutations
      Mut_codpos_1 = nrow(MEGA_CDS[MEGA_CDS$Codon_Pos == 1,])
      Mut_codpos_2 = nrow(MEGA_CDS[MEGA_CDS$Codon_Pos == 2,])
      Mut_codpos_3 = nrow(MEGA_CDS[MEGA_CDS$Codon_Pos == 3,])
      ## Computing proportions
      Mut_prop_cospos_1 = Mut_codpos_1 / sum(Mut_codpos_1,Mut_codpos_2,Mut_codpos_3)
      Mut_prop_cospos_2 = Mut_codpos_2 / sum(Mut_codpos_1,Mut_codpos_2,Mut_codpos_3)
      Mut_prop_cospos_3 = Mut_codpos_3 / sum(Mut_codpos_1,Mut_codpos_2,Mut_codpos_3)
    }
    
    ##### Per Site Degeneracy #####
    {
      ## Counting number of mutations
      Mut_0x = nrow(MEGA_0x)
      Mut_2x = nrow(MEGA_2x)
      Mut_4x = nrow(MEGA_4x)
      ## Computing proportions
      # here I need to count the number of 0x, 2x and 4x sites
      # Load the REF genome to count on, only CDS sequence
      REF_CDS = read.fasta("../01-DATA/0-REF/REF_CDS.fasta",forceDNAtolower = F)
      # Use the function to tag each site with their degeneracy
      degeneracy_tag_CDS(REF_CDS)
      # Count each in a table
      Degeneracy_per_site <- data.frame(matrix(ncol=4))
      colnames(Degeneracy_per_site) <- c("Gene","Zerofold","Twofold","Fourfold")
      for (gene in names(Seq_table_list)) {
        gene_table <- data.frame(Seq_table_list[gene])
        colnames(gene_table) <- c("Seq","Pos","Degeneracy")
        Zero <- nrow(gene_table[gene_table$Degeneracy == "Zerofold",])
        Two <- nrow(gene_table[gene_table$Degeneracy == "Twofold",])
        Four <- nrow(gene_table[gene_table$Degeneracy == "Fourfold",])
        Degeneracy_per_site <- rbind(Degeneracy_per_site,c(gene,Zero,Two,Four))
      }
      Degeneracy_per_site <- Degeneracy_per_site[2:13,]
      # Write this table to use again if necessary
      write.table(Degeneracy_per_site,file = "R_txt_files/Degeneracy_sites_per_gene.txt",sep="\t",quote = F,row.names = F)
      # Finally computing the proportions
      Mut_0x_prop <- Mut_0x / sum(as.numeric(Degeneracy_per_site$Zerofold))
      Mut_2x_prop <- Mut_2x / sum(as.numeric(Degeneracy_per_site$Twofold))
      Mut_4x_prop <- Mut_4x / sum(as.numeric(Degeneracy_per_site$Fourfold))
    }
  }
  
  ##### Ts/Tv #####
  {
    Transi <- MEGA[MEGA$Transition == "Transition",]
    Transver <- MEGA[MEGA$Transition == "Transversion",]
    ##### Per Codon position #####
    {
      Tsall <- nrow(Transi)
      Ts1 <- nrow(Transi[Transi$Codon_Pos == 1,])
      Ts2 <- nrow(Transi[Transi$Codon_Pos == 2,])
      Ts3 <- nrow(Transi[Transi$Codon_Pos == 3,])
      
      Tvall <- nrow(Transver)
      Tv1 <- nrow(Transver[Transver$Codon_Pos == 1,])
      Tv2 <- nrow(Transver[Transver$Codon_Pos == 2,])
      Tv3 <- nrow(Transver[Transver$Codon_Pos == 3,])
      
      TsTvall = Tsall/Tvall
      TsTv1 = Ts1/Tv1
      TsTv2 = Ts2/Tv2
      TsTv3 = Ts3/Tv3
    }
    
    ##### Per Site Degeneracy #####
    {
      Ts0x <- nrow(Transi[Transi$Site_4x_2x == "0x",])
      Ts2x <- nrow(Transi[Transi$Site_4x_2x == "2x",])
      Ts4x <- nrow(Transi[Transi$Site_4x_2x == "4x",])
      
      Tv0x <- nrow(Transver[Transver$Site_4x_2x == "0x",])
      Tv2x <- nrow(Transver[Transver$Site_4x_2x == "2x",])
      Tv4x <- nrow(Transver[Transver$Site_4x_2x == "4x",])
      
      TsTv0x <- Ts0x/Tv0x
      TsTv2x <- Ts2x/Tv2x
      TsTv4x <- Ts4x/Tv4x
      
      # Intergenic
      TsInter <- nrow(MEGA_NON_CODING[MEGA_NON_CODING$Transition == "Transversion",])
      TvInter <- nrow(MEGA_NON_CODING[MEGA_NON_CODING$Transition == "Transition",])
      TsTvInter <- TsInter/TvInter
    }
  }
}

##### Statistical tests #####
{
  # Make vector of mut counts for each coding position
  mut_codpos <- c(Mut_codpos_1,Mut_codpos_2,Mut_codpos_3)
  
  # Waneka data
  Waneka = read.table("Previous_results/Waneka_my_annotation.txt",header = T,sep="\t")
  Waneka_CDS <- Waneka[Waneka$Gene_Type == "CDS",]
  table(Waneka_CDS$Codon_position) / sum(table(Waneka_CDS$Codon_position))
  Waneka_mut_codpos <- c(70,86,49) 
  
  # Expected is neutral
  expected <- c(1/3,1/3,1/3)
  
  # G.tests of the counts vs expected proportions under neutrality
  g.test(mut_codpos,p=expected)
  g.test(Waneka_mut_codpos,p=expected)
}