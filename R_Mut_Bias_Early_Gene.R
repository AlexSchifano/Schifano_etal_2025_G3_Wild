setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")

## I want to look into the number of mutations in the early part of genes
## Also check AT content, syn-site proportion
## Ultimately I have to identify if there is a bias for less syn mutations in early gene and if there is and explanation


########### 1) Mutation per codon ###########
{
  # Load CDS mutations
  MEGA_CDS <- read.table("R_txt_files/Mutation_CDS.txt", header = T, sep = "\t")
  MEGA_CDS <- MEGA_CDS[MEGA_CDS$Synonymous == "Synonymous",]
  
  # Load the GTF to get gene start info and subset CDS
  GTF <- read.table("../01-DATA/0-REF/WS283_MtDNA.gtf", col.names = c("Chr","START","END","Gene_ID","Gene_type"))
  GTF <- GTF[GTF$Gene_type == "CDS",]
  GTF$SIZE <- GTF$END - GTF$START + 1
  
  ## I use the pos ref to compute the number of mutation per codon position ==> pos 1-3 = 1 ; pos 4-6 = 2 etc
  ## I store everything in this table
  # Cox1 being the largest gene I use it for codon count
  nb_codon_cox1 <- seq(from = 1, to = GTF$SIZE[GTF$Gene_ID == "COX1"] / 3)
  Mut_per_Codon <- data.frame(Codon_nb = nb_codon_cox1, Number_gene = 0, Mutation_Counts = 0)
  # Adding number of genes with a codon at this codon position
  for (ge in GTF$Gene_ID) {
    nb_codon <- seq(from = 1, to = GTF$SIZE[GTF$Gene_ID == ge] / 3)
    Mut_per_Codon$Number_gene[Mut_per_Codon$Codon_nb %in% nb_codon] <- Mut_per_Codon$Number_gene[Mut_per_Codon$Codon_nb %in% nb_codon] + 1
  }
  
  ## Loop to count the mutations per codon position
  for (ro in 1:nrow(MEGA_CDS)) {
    pos_ref <- MEGA_CDS[ro,"POS_REF"]
    gene <- MEGA_CDS[ro,"Gene_ID"]
    pos_gene = pos_ref - GTF$START[GTF$Gene_ID == gene] + 1
    pos_codon = ceiling(pos_gene / 3)
    Mut_per_Codon$Mutation_Counts[Mut_per_Codon$Codon_nb == pos_codon] <- Mut_per_Codon$Mutation_Counts[Mut_per_Codon$Codon_nb == pos_codon] + 1
    MEGA_CDS[ro,"POS_Codon"] <- pos_codon
    MEGA_CDS[ro,"POS_Gene"] <- pos_gene
    }
  
  ## Pulling genes together
  {
    ## Plot all
    library(ggplot2)
    ggplot(Mut_per_Codon, aes(x=Codon_nb,y=Mutation_Counts)) + geom_line() + geom_smooth()
    
    ## Plot subset first 78 positions (length of shortest)
    Mut_per_Codon_78 <- Mut_per_Codon[Mut_per_Codon$Codon_nb < 79,]
    ggplot(Mut_per_Codon_78, aes(x=Codon_nb,y=Mutation_Counts)) + geom_line() + geom_smooth()
    
    ## Same but grouping codons by group of 3
    # Example of grouping codons into 9 positions
    Mut_Codon_triplet <- Mut_per_Codon_78
    Mut_Codon_triplet$Triplet_group <- rep(1:ceiling(nrow(Mut_per_Codon_78)/3), each = 3)[1:nrow(Mut_per_Codon_78)]
    library(dplyr)
    
    # Sum mutations per triplet group
    Triplet_counts <- Mut_Codon_triplet %>%
      group_by(Triplet_group) %>%
      summarise(Triplet_mutations = sum(Mutation_Counts))
    
    # Subset 20 * 9 nucleotides = 180bp
    Triplet_counts <- Triplet_counts[Triplet_counts$Triplet_group < 21,]
    # Plot
    Mut <- ggplot(Triplet_counts, aes(x=Triplet_group,y=Triplet_mutations)) + geom_line() + geom_smooth() + theme_minimal()
    Mut
    
    ggsave(Mut, filename = "../03-SV_ANALYSES/Neighbour_Effect/Early_gene_Bias/RAW_Fig/Fig_Mut_Early.pdf", height = 15, width = 15)
  }
  
  ## NOT WORKING GAVE UP ## With count per gene 
  {
    # Table of gene
    Gene_codon_count <- data.frame(Pos = 1:(max(GTF$SIZE)), ND6=0, ND4L=0, ND1=0,
                                   ATP6=0,ND2=0,CTB1=0,COX3=0,ND4=0,COX1=0,COX2=0,ND3=0,ND5=0)
    # Create groups per position
    Gene_codon_count$Codon_Pos <- ceiling(Gene_codon_count$Pos / 3)
    Gene_codon_count$Group_9 <- ceiling(Gene_codon_count$Pos / 9)
    
    # Loop to count mut per pos
    for (ro in 1:nrow(MEGA_CDS)){
      gene <- MEGA_CDS[ro,"Gene_ID"]
      pos_gene <- MEGA_CDS[ro,"POS_Gene"]
      Gene_codon_count[pos_gene,gene] <- Gene_codon_count[pos_gene,gene] + 1
    }
    
    # Group by group of 9 nuc
    {
      Gene_count_9 <- data.frame(Group_9 = unique(Gene_codon_count$Group_9), ND6=0, ND4L=0, ND1=0,
                                 ATP6=0,ND2=0,CTB1=0,COX3=0,ND4=0,COX1=0,COX2=0,ND3=0,ND5=0)
      gene_list <- c("ND6","ND4L","ND1","ATP6","ND2","CTB1","COX3","ND4","COX1","COX2","ND3","ND5")
      for (grp in Gene_count_9$Group_9){
        Gene_count_9[grp,"ND6"] <- sum(Gene_codon_count$ND6[Gene_codon_count$Group_9 == grp])
        Gene_count_9[grp,"ND4L"] <- sum(Gene_codon_count$ND4L[Gene_codon_count$Group_9 == grp])
        Gene_count_9[grp,"ND1"] <- sum(Gene_codon_count$ND1[Gene_codon_count$Group_9 == grp])
        Gene_count_9[grp,"ATP6"] <- sum(Gene_codon_count$ATP6[Gene_codon_count$Group_9 == grp])
        Gene_count_9[grp,"ND2"] <- sum(Gene_codon_count$ND2[Gene_codon_count$Group_9 == grp])
        Gene_count_9[grp,"CTB1"] <- sum(Gene_codon_count$CTB1[Gene_codon_count$Group_9 == grp])
        Gene_count_9[grp,"COX3"] <- sum(Gene_codon_count$COX3[Gene_codon_count$Group_9 == grp])
        Gene_count_9[grp,"ND4"] <- sum(Gene_codon_count$ND4[Gene_codon_count$Group_9 == grp])
        Gene_count_9[grp,"COX1"] <- sum(Gene_codon_count$COX1[Gene_codon_count$Group_9 == grp])
        Gene_count_9[grp,"COX2"] <- sum(Gene_codon_count$COX2[Gene_codon_count$Group_9 == grp])
        Gene_count_9[grp,"ND3"] <- sum(Gene_codon_count$ND3[Gene_codon_count$Group_9 == grp])
        Gene_count_9[grp,"ND5"] <- sum(Gene_codon_count$ND5[Gene_codon_count$Group_9 == grp])
      }
    }
    
    # Plot
    Gene_count_9 <- Gene_count_9[1:20,]
    
    ggplot() + geom_point(data=Gene_count_9, aes(x=Group_9 ,y=ND6)) + geom_point(data=Gene_count_9, aes(x=Group_9 ,y=ND4L)) +
      geom_point(data=Gene_count_9, aes(x=Group_9 ,y=ND1)) + geom_point(data=Gene_count_9, aes(x=Group_9 ,y=ATP6)) +
      geom_point(data=Gene_count_9, aes(x=Group_9 ,y=ND2)) + geom_point(data=Gene_count_9, aes(x=Group_9 ,y=CTB1)) +
      geom_point(data=Gene_count_9, aes(x=Group_9 ,y=COX3)) + geom_point(data=Gene_count_9, aes(x=Group_9 ,y=ND4)) +
      geom_point(data=Gene_count_9, aes(x=Group_9 ,y=COX1)) + geom_point(data=Gene_count_9, aes(x=Group_9 ,y=COX2)) +
      geom_point(data=Gene_count_9, aes(x=Group_9 ,y=ND3)) + geom_point(data=Gene_count_9, aes(x=Group_9 ,y=ND5))
    

    
  }
}

########### 2) AT-content ###########
{
  library(seqinr)
  WS283_CDS <- read.fasta("../01-DATA/0-REF/FULL_REF/WS283_MtDNA_CDS.fasta", forceDNAtolower = F)
  
  ## I focus on the first 78 codon = 234 bp (ND4L size)
  ## I extract each sequence in a table
  AT_Count <- data.frame(Pos = seq(1:234),ND6="-",ND4L="-",ND1="-",ATP6="-",ND2="-",CTB1="-",COX3="-",ND4="-",COX1="-",COX2="-",ND3="-",ND5="-")
  # Gene order should be the same between this and the fasta file
  for (gen in 1:12) {
    AT_Count[1:234,gen+1] <- WS283_CDS[[gen]][1:234]
  }
  
  ## I compute AT % per position
  AT_Count$AT_perc <- "-"
  for (ro in seq(from = 1, to = nrow(AT_Count))){
    pos <- AT_Count[ro,2:13]
    AT_perc <- (sum(pos %in% c("A","T")) / 12) * 100
    AT_Count[ro,"AT_perc"] <- as.numeric(AT_perc)
  }
  
  ## Plot AT per pos
  ggplot(AT_Count, aes(x=as.numeric(Pos), y=as.numeric(AT_perc))) + geom_line() + geom_smooth()
  
  # Same but doing similar grouping as triplet of codon ==> per 9bp
  AT_count_nine <- AT_Count
  AT_count_nine$Group <- rep(1:ceiling(nrow(AT_Count)/9), each = 9)[1:nrow(AT_Count)]
  library(dplyr)
  
  # Sum mutations per triplet group
  AT_count_nine$AT_perc <- as.numeric(AT_count_nine$AT_perc)
  AT_nine_counts <- AT_count_nine %>%
    group_by(Group) %>%
    summarise(Mean_AT_perc = mean(AT_perc, na.rm = TRUE))  # Add na.rm = TRUE to handle any missing values
  
  # Subset first 20*9 nuc
  AT_nine_counts <- AT_nine_counts[AT_nine_counts$Group < 21,]
  
  ## Plot AT per pos
  AT <- ggplot(AT_nine_counts, aes(x = as.numeric(Group), y = as.numeric(Mean_AT_perc))) +
    geom_line() +
    geom_smooth() +
    theme_minimal() +
    coord_cartesian(ylim = c(65, 90))  # Adjust based on actual data range
  AT
  
  ggsave(AT, filename = "../03-SV_ANALYSES/Neighbour_Effect/Early_gene_Bias/RAW_Fig/Fig_AT_Early.pdf", height = 15, width = 15)
}

########### 3) Mutation per codon with distance to STOP ###########
{
  # Load CDS mutations
  MEGA_CDS <- read.table("R_txt_files/Mutation_CDS.txt", header = T, sep = "\t")
  MEGA_CDS <- MEGA_CDS[MEGA_CDS$Synonymous == "Synonymous",]
  
  # Load the GTF to get gene start info and subset CDS
  GTF <- read.table("../01-DATA/0-REF/WS283_MtDNA.gtf", col.names = c("Chr","START","END","Gene_ID","Gene_type"))
  GTF <- GTF[GTF$Gene_type == "CDS",]
  GTF$SIZE <- GTF$END - GTF$START + 1
  
  ## I use the pos ref to compute the number of mutation per codon position ==> pos 1-3 = 1 ; pos 4-6 = 2 etc
  ## I store everything in this table
  Mut_dist_stop <- data.frame(Dist_stop = seq(1:max(GTF$SIZE)), Number_gene = 12, Mutation_Counts = 0)
  
  ## Loop to count the mutations per codon position
  for (ro in 1:nrow(MEGA_CDS)) {
    pos_ref <- MEGA_CDS[ro,"POS_REF"]
    gene <- MEGA_CDS[ro,"Gene_ID"]
    dist_stop <- GTF$END[GTF$Gene_ID == gene] - pos_ref
    Mut_dist_stop$Mutation_Counts[Mut_dist_stop$Dist_stop == dist_stop] <- Mut_dist_stop$Mutation_Counts[Mut_dist_stop$Dist_stop == dist_stop] + 1
  }
  
  ## Same but grouping nucleotides by groups of 9
  # Example of grouping 
  # Subset last 200bp (208 because 9*23 = 207)
  Mut_dist_stop_nine <- Mut_dist_stop[Mut_dist_stop$Dist_stop < 208,]
  Mut_dist_stop_nine$Group <- rep(1:ceiling(nrow(Mut_dist_stop_nine)/9), each = 9)[1:nrow(Mut_dist_stop_nine)]
  
  # Sum mutations per triplet group
  library(dplyr)
  Group_counts <- Mut_dist_stop_nine %>%
    group_by(Group) %>%
    summarise(Triplet_mutations = sum(Mutation_Counts))
  
  # Subset first 20 groups of 9
  Group_counts <- Group_counts[Group_counts$Group < 21,]
  
  # Plot
  Mut_Stop <- ggplot(Group_counts, aes(x=Group,y=as.numeric(Triplet_mutations))) + geom_line() + geom_smooth() + theme_minimal() +
    xlab("Distance from stop x9") + ylab("Mutation count per 9 nuc")
  Mut_Stop
  
  ggsave(Mut_Stop, filename = "../03-SV_ANALYSES/Neighbour_Effect/Early_gene_Bias/RAW_Fig/Fig_Mut_Late.pdf", height = 15, width = 15)

}

########### 4) AT-content close to stop ###########
{
  library(seqinr)
  WS283_CDS <- read.fasta("../01-DATA/0-REF/FULL_REF/WS283_MtDNA_CDS.fasta", forceDNAtolower = F)
  
  ## I focus on the first 78 codon = 234 bp (ND4L size)
  ## I extract each sequence in a table
  AT_Count <- data.frame(Dist_stop = seq(1:207),ND6="-",ND4L="-",ND1="-",ATP6="-",ND2="-",CTB1="-",COX3="-",ND4="-",COX1="-",COX2="-",ND3="-",ND5="-")
  # Gene order should be the same between this and the fasta file
  for (gen in 1:12) {
    AT_Count[1:207,gen+1] <- WS283_CDS[[gen]][seq(from = length(WS283_CDS[[gen]]), to = length(WS283_CDS[[gen]]) - 206)]
  }
  
  ## I compute AT % per position
  AT_Count$AT_perc <- "-"
  for (ro in seq(from = 1, to = nrow(AT_Count))){
    pos <- AT_Count[ro,2:13]
    AT_perc <- (sum(pos %in% c("A","T")) / 12) * 100
    AT_Count[ro,"AT_perc"] <- as.numeric(AT_perc)
  }
  
  ## Plot AT per pos
  ggplot(AT_Count, aes(x=as.numeric(Dist_stop), y=as.numeric(AT_perc))) + geom_line() + geom_smooth()
  
  # Same but doing similar grouping as triplet of codon ==> per 9bp
  AT_count_nine <- AT_Count
  AT_count_nine$Group <- rep(1:ceiling(nrow(AT_Count)/9), each = 9)[1:nrow(AT_Count)]
  library(dplyr)
  
  # Sum mutations per triplet group
  AT_count_nine$AT_perc <- as.numeric(AT_count_nine$AT_perc)
  AT_nine_counts <- AT_count_nine %>%
    group_by(Group) %>%
    summarise(Mean_AT_perc = mean(AT_perc, na.rm = TRUE))  # Add na.rm = TRUE to handle any missing values
  
  # Subset fist 20 groups
  AT_nine_counts <- AT_nine_counts[AT_nine_counts$Group < 21,]
  ## Plot AT per pos
  AT_STOP <- ggplot(AT_nine_counts, aes(x = as.numeric(Group), y = as.numeric(Mean_AT_perc))) +
    geom_line() + xlab("Distance from STOP") + ylab("AT-content") +
    geom_smooth() +
    theme_minimal() +
    coord_cartesian(ylim = c(65, 90))  # Adjust based on actual data range
  
  AT_STOP
  
  ggsave(AT_STOP, filename = "../03-SV_ANALYSES/Neighbour_Effect/Early_gene_Bias/RAW_Fig/Fig_AT_Late.pdf", height = 15, width = 15)
}

########## 5) AT-content thoughout the mtDNA ##########
{

  # Function to look at a specific region
  AT_Perc_Region <- function(START = 0, END, Global_Avg =T) {
    ## Compute Genome wide AT-content
    {
      # Read ref
      REF <- read.fasta("~/Documents/C.elegans_sync/Final_files/01-DATA/0-REF/WS283_MtDNA.fasta",forceDNAtolower = F)
      # put as dataframe
      REF_df <- data.frame(pos = 1:length(REF[[1]]), Seq = REF[[1]])
      # Add groups of 9 nucleotides
      REF_df$Group <- ceiling(REF_df$pos / 9)
      # Group AT-content
      Group_AT_REF <- data.frame(Group = unique(REF_df$Group), AT_Perc = 0)
      for (grp in unique(REF_df$Group)){
        AT_perc <- (sum(REF_df$Seq[REF_df$Group == grp] == "A") + sum(REF_df$Seq[REF_df$Group == grp] == "T")) / 9
        Group_AT_REF$AT_Perc[Group_AT_REF$Group == grp] <- AT_perc
      }
      # Compute global average
      Avg_AT <- mean(Group_AT_REF$AT_Perc)
      sd_AT <- sd(Group_AT_REF$AT_Perc)
    }
    
    ## Subset the table based on START and END
    {
      # Set END
      if (exists("END") == F) { END = length(REF[[1]]) }
      
      # Set it back to match the 9 nucleotide window size
      START_9 <- ceiling(START/9)
      END_9 <- ceiling(END/9)
      
      # Subset the table
      Group_AT_REF_Sub <- Group_AT_REF[Group_AT_REF$Group %in% seq(from = START_9, to = END_9),]
      
      # Write to global env
      name <- paste("AT_Perc_9bp_",START,"_",END,sep="")
      assign(name,Group_AT_REF_Sub, envir = globalenv())
      
      # Plot
      if (Global_Avg) {
        Plot <- ggplot(Group_AT_REF_Sub, aes(x=Group,y=AT_Perc)) + geom_line() + 
          geom_smooth() + theme_minimal() + xlab("Position x9") + 
          geom_hline(yintercept = Avg_AT, col="red") +
          geom_hline(yintercept = Avg_AT + sd_AT, col = "red", linetype="dashed") +
          geom_hline(yintercept = Avg_AT - sd_AT, col = "red", linetype="dashed")
      } else {
        Plot <- ggplot(Group_AT_REF_Sub, aes(x=Group,y=AT_Perc)) + geom_line() + 
          geom_smooth() + theme_minimal() + xlab("Position x9")
      }

      
      print(Plot)
    
    }
  }
  
  # Test
  AT_Perc_Region(START = 113, END = 300, Global_Avg = F)
  
  # Before ctb-1 del (mtDNA:4809-5307)
  AT_Perc_Region(START = 4700, END = 4900)
  # After ctb-1 del
  AT_Perc_Region(START = 5200, END = 5400)
  # In the deletion region
  AT_Perc_Region(START = 4809, END = 5307)
  
  # Before cox3/nd4 del (mtDNA:6361-7394)
  AT_Perc_Region(START = 6000, END = 6361)
  # After cox3/nd4  del
  AT_Perc_Region(START = 7394, END = 7700)
  # In the deletion region
  AT_Perc_Region(START = 7300, END = 7450)
  
  # Before cox1 del (mtDNA:7920-9104)
  AT_Perc_Region(START = 7600, END = 7920)
  # After cox1  del
  AT_Perc_Region(START = 9104, END = 9400)
  # spanning breakpoints
  AT_Perc_Region(START = 7800, END = 8000)
  
  
}

########## 6) Correlations Syn X GC ##########
{
  #### Syn Mut X GC_1-2 content ####
  {
    ### GC-content
    {
      ## GC_1-2 content count along genome
      library(seqinr)
      WS283_CDS <- read.fasta("../01-DATA/0-REF/FULL_REF/WS283_MtDNA_CDS.fasta", forceDNAtolower = F)
      
      ## I focus on the first 78 codon = 234 bp (ND4L size)
      ## I extract each sequence in a table
      GC_Count <- data.frame(Pos = seq(1:1584),ND6="-",ND4L="-",ND1="-",ATP6="-",ND2="-",CTB1="-",COX3="-",ND4="-",COX1="-",COX2="-",ND3="-",ND5="-")
      for (gen in 1:12) {
        GC_Count[1:1584,gen+1] <- WS283_CDS[[gen]][1:1584]
      }
      
      ## Add nb gene with pos
      for (pos in GC_Count$Pos){
        GC_Count[pos,"Nb_Gene"] <- 12 - sum(is.na(GC_Count[pos,2:13]))
      }
      
      ## I compute GC % per position
      GC_Count$GC_perc <- "-"
      for (ro in seq(from = 1, to = nrow(GC_Count))){
        pos <- GC_Count[ro,2:13]
        pos <- pos[is.na(pos) == F]
        GC_perc <- 1 - ((sum(pos %in% c("A","T")) / length(pos)))
        GC_Count[ro,"GC_perc"] <- as.numeric(GC_perc)
      }
      
      ## Add Codon_pos
      GC_Count$Codon_pos <- rep(c(1,2,3))
      
      ## Add group for later
      GC_Count$Group <- ceiling(GC_Count$Pos / 3)
    }
    
    ### Mut per gene pos
    {
      # Load CDS mutations
      MEGA_CDS <- read.table("R_txt_files/Mutation_CDS.txt", header = T, sep = "\t")
      MEGA_CDS <- MEGA_CDS[MEGA_CDS$Synonymous == "Synonymous",]
      
      # Load the GTF to get gene start info and subset CDS
      GTF <- read.table("../01-DATA/0-REF/WS283_MtDNA.gtf", col.names = c("Chr","START","END","Gene_ID","Gene_type"))
      GTF <- GTF[GTF$Gene_type == "CDS",]
      GTF$SIZE <- GTF$END - GTF$START + 1
      
      ## I use the pos ref to compute the number of mutation per codon position ==> pos 1-3 = 1 ; pos 4-6 = 2 etc
      ## I store everything in this table
      # Cox1 being the largest gene I use it for codon count
      nb_codon_cox1 <- seq(from = 1, to = GTF$SIZE[GTF$Gene_ID == "COX1"] / 3)
      Mut_per_Codon <- data.frame(Codon_nb = nb_codon_cox1, Number_gene = 0, Mutation_Counts = 0)
      # Adding number of genes with a codon at this codon position
      for (ge in GTF$Gene_ID) {
        nb_codon <- seq(from = 1, to = GTF$SIZE[GTF$Gene_ID == ge] / 3)
        Mut_per_Codon$Number_gene[Mut_per_Codon$Codon_nb %in% nb_codon] <- Mut_per_Codon$Number_gene[Mut_per_Codon$Codon_nb %in% nb_codon] + 1
      }
      
      ## Loop to count the mutations per codon position
      for (ro in 1:nrow(MEGA_CDS)) {
        pos_ref <- MEGA_CDS[ro,"POS_REF"]
        gene <- MEGA_CDS[ro,"Gene_ID"]
        pos_gene = pos_ref - GTF$START[GTF$Gene_ID == gene] + 1
        pos_codon = ceiling(pos_gene / 3)
        Mut_per_Codon$Mutation_Counts[Mut_per_Codon$Codon_nb == pos_codon] <- Mut_per_Codon$Mutation_Counts[Mut_per_Codon$Codon_nb == pos_codon] + 1
        MEGA_CDS[ro,"POS_Codon"] <- pos_codon
        MEGA_CDS[ro,"POS_Gene"] <- pos_gene
      }
    }
    
    ### For 2x sites ###
    {
      MEGA_2x <- MEGA_CDS[MEGA_CDS$Site_4x_2x == "2x",]
      
      ### Add the mut counts to the GC table
      for (pos in GC_Count$Pos){
        GC_Count[pos,"Mut_Count"] <- nrow(MEGA_2x[MEGA_2x$POS_Gene == pos,])
      }
      
      # Sum mutations per group 
      Mut_group <- GC_Count %>%
        group_by(Group) %>%
        summarise(Mut_Count = sum(Mut_Count), Nb_Gene = mean(Nb_Gene)) 
      
      ## Subset position 1 and 2
      GC_12_Count <- GC_Count[GC_Count$Codon_pos != 3,]
      # GC-content per group
      GC_12_Group <- GC_12_Count %>%
        group_by(Group) %>%
        summarise(GC_12 = mean(as.numeric(GC_perc)))
      # Add the mut counts
      GC_12_Group_Mut <- cbind(GC_12_Group,Mut_group[,-1])
      GC_12_Group_Mut$GC_12 <- round(GC_12_Group_Mut$GC_12,3)
      # Compute Pearson correlation
      cor_test <- cor.test(GC_12_Group_Mut$GC_12, GC_12_Group_Mut$Mut_Count)
      r_value <- round(cor_test$estimate, 4)
      p_value <- signif(cor_test$p.value, 4)
      # Plot
      ggplot(GC_12_Group_Mut, aes(x = factor(GC_12), y = Mut_Count)) +
        geom_boxplot() +
        geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "blue") +
        theme_minimal() + ggtitle("Two-fold mutation X GC12") + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
        ylab("Mutation Count") + xlab("GC12") 
      
      
      ## Subset pos 3
      GC_3_Count <- GC_Count[GC_Count$Codon_pos == 3,]
      # GC-content per group (group = codon)
      GC_3_Group <- GC_3_Count %>%
        group_by(Group) %>%
        summarise(GC_3 = mean(as.numeric(GC_perc)))
      # Add the mut counts
      GC_3_Group_Mut <- cbind(GC_3_Group,Mut_group[,-1])
      GC_3_Group_Mut$GC_3 <- round(GC_3_Group_Mut$GC_3,3)
      # Plot
      ggplot(GC_3_Group_Mut, aes(x = factor(GC_3), y = Mut_Count)) +
        geom_boxplot() +
        geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "blue") +
        theme_minimal() + ggtitle("Two-fold mutation X GC3") + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
        ylab("Mutation Count") + xlab("GC3")
      # Correlation
      cor(GC_3_Group_Mut$GC_3,GC_3_Group_Mut$Mut_Count, method = "pearson") 
      cor.test(GC_3_Group_Mut$GC_3,GC_3_Group_Mut$Mut_Count, method = "pearson")
      mod_2X_GC3 <- lm(GC_3_Group_Mut$GC_3 ~ GC_3_Group_Mut$Mut_Count)
      summary(mod_2X_GC3)
    }
    
    ### For 4x sites ###
    {
      MEGA_4x <- MEGA_CDS[MEGA_CDS$Site_4x_2x == "4x",]
      
      ### Add the mut counts to the GC table
      for (pos in GC_Count$Pos){
        GC_Count[pos,"Mut_Count"] <- nrow(MEGA_4x[MEGA_4x$POS_Gene == pos,])
      }
      
      # Sum mutations per group 
      Mut_group <- GC_Count %>%
        group_by(Group) %>%
        summarise(Mut_Count = sum(Mut_Count), Nb_Gene = mean(Nb_Gene)) 
      
      ## Subset position 1 and 2
      GC_12_Count <- GC_Count[GC_Count$Codon_pos != 3,]
      # GC-content per group
      GC_12_Group <- GC_12_Count %>%
        group_by(Group) %>%
        summarise(GC_12 = mean(as.numeric(GC_perc)))
      # Add the mut counts
      GC_12_Group_Mut <- cbind(GC_12_Group,Mut_group[,-1])
      GC_12_Group_Mut$GC_12 <- round(GC_12_Group_Mut$GC_12,digits = 3)
      # Plot
      ggplot(GC_12_Group_Mut, aes(x = factor(GC_12), y = Mut_Count)) +
        geom_boxplot() +
        geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "blue") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
        ggtitle("Four-fold mutations * GC12") + ylab("Mutation Count") + xlab("GC12")
      # Correlation
      cor(GC_12_Group_Mut$GC_12,GC_12_Group_Mut$Mut_Count, method = "pearson") 
      cor.test(GC_12_Group_Mut$GC_12,GC_12_Group_Mut$Mut_Count, method = "pearson")
      mod_4x_GC12 <- lm(GC_12_Group_Mut$GC_12 ~ GC_12_Group_Mut$Mut_Count)
      summary(mod_4x_GC12)
      library(olsrr)
      ols_plot_resid_qq(mod_4x_GC12)
      ols_plot_resid_fit(mod_4x_GC12)
      
      ## Subset pos 3
      GC_3_Count <- GC_Count[GC_Count$Codon_pos == 3,]
      # GC-content per group (group = codon)
      GC_3_Group <- GC_3_Count %>%
        group_by(Group) %>%
        summarise(GC_3 = mean(as.numeric(GC_perc)))
      # Add the mut counts
      GC_3_Group_Mut <- cbind(GC_3_Group,Mut_group[,-1])
      GC_3_Group_Mut$GC_3 <- round(GC_3_Group_Mut$GC_3,3)
      # Plot
      ggplot(GC_3_Group_Mut, aes(x = factor(GC_3), y = Mut_Count)) +
        geom_boxplot() +
        geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "blue") +
        theme_minimal() + ggtitle("Four-fold mutation X GC3") + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
        ylab("Mutation Count") + xlab("GC3")
      # Correlation
      cor(GC_3_Group_Mut$GC_3,GC_3_Group_Mut$Mut_Count, method = "pearson") 
      cor.test(GC_3_Group_Mut$GC_3,GC_3_Group_Mut$Mut_Count, method = "pearson")
      mod_4x_GC3 <- lm(GC_3_Group_Mut$GC_3 ~ GC_3_Group_Mut$Mut_Count)
      summary(mod_4x_GC3)
      library(olsrr)
      ols_plot_resid_qq(mod_4x_GC3)
      ols_plot_resid_fit(mod_4x_GC3)
    }
    
    ### Do GC12 x Nb_Mut per gene (12 data points)
    {
      ## GC per gene
      {
        # Subsets
        GC12_count <- GC_Count[GC_Count$Codon_pos %in% c(1,2),]
        GC3_count <- GC_Count[GC_Count$Codon_pos == 3,]
        MEGA <- read.table("R_txt_files/Mutation_List.txt", header = T, sep="\t")
        # Setting up the table
        GC_gene <- data.frame(row.names = colnames(GC_Count[1,2:13]), GC12_perc=rep(0,12), GC3_perc=rep(0,12), Fourfold_mut=rep(0,12), Twofold_mut=rep(0,12), Syn_mut=rep(0,12))
        # Loop for counts per gene
        for (g in row.names(GC12_gene)) {
          GC_gene[g,"GC12_perc"] <- sum(GC12_count[,g] %in% c("G","C")) / sum(GC12_count[,g] %in% c("G","C","A","T"))
          GC_gene[g,"GC3_perc"] <- sum(GC3_count[,g] %in% c("G","C")) / sum(GC3_count[,g] %in% c("G","C","A","T"))
          MEGA_g <- MEGA[MEGA$Gene_ID == g,]
          GC_gene[g,"Fourfold_mut"] <- nrow(MEGA_g[MEGA_g$Site_4x_2x == "4x",]) 
          GC_gene[g,"Twofold_mut"] <- nrow(MEGA_g[MEGA_g$Site_4x_2x == "2x",]) 
          GC_gene[g,"Syn_mut"] <- nrow(MEGA_g[MEGA_g$Synonymous == "Synonymous",]) 
          GC_gene[g,"Size"] <- GTF$SIZE[GTF$Gene_ID == g]
        }
        # Adding number of Syn sites for normalization as well as 4x sites
        {
          # In the same gene order
          GC_gene$Syn_Sites <- c(91,43.3333,185.3333,127.3333,167.66667,234,160.33333,
                                 257.33333,357.66667,142.3333,64.66667,326)
          
          Degeneracy_gene <- read.table("R_txt_files/Degeneracy_sites_per_gene.txt",header = T,sep="\t")
          GC_gene$Fourfold_sites <- Degeneracy_gene$Fourfold
          # Normalized
          GC_gene$Normalised_Syn <- GC_gene$Syn_mut / GC_gene$Syn_Sites
          GC_gene$Normalised_Fourfold <- GC_gene$Fourfold_mut / GC_gene$Fourfold_sites
        }

      }
      
      ## Testing correlation
      {
        ## Normalised_Syn
        {
          # GC12 X Normalised_Syn
          ggplot(GC_gene, aes(x=GC12_perc,y=Normalised_Syn)) + geom_point() + geom_smooth(method = "lm") + theme_minimal() +
            ylab("Normalised Synonymous mutations") + xlab("GC12")
          # Test correlation
          cor.test(GC_gene$GC12_perc,GC_gene$Normalised_Syn)
          mod_GC12_Syn <- lm(data = GC_gene, GC12_perc ~ Normalised_Syn)
          summary(mod_GC12_Syn)
          
          # GC3 X Normalised_Syn
          ggplot(GC_gene, aes(x=GC3_perc,y=Normalised_Syn)) + geom_point() + geom_smooth(method = "lm") + theme_minimal() +
            ylab("Number of Synonymous mutations") + xlab("GC3")
          # Test correlation
          cor.test(GC_gene$GC3_perc,GC_gene$Normalised_Syn)
          mod_GC3_Syn <- lm(data = GC_gene, GC3_perc ~ Normalised_Syn)
          summary(mod_GC3_Syn)
        }
        
        ## Normalised Fourfold mut
        {
          # GC12 X 4x_mut
          ggplot(GC_gene, aes(x=GC12_perc,y=Normalised_Fourfold)) + geom_point() + geom_smooth(method = "lm") + theme_minimal() +
            ylab("Normalised Fourfold mutations") + xlab("GC12")
          # Test correlation
          cor.test(GC_gene$GC12_perc,GC_gene$Normalised_Fourfold)
          mod_GC12_4x <- lm(data = GC_gene, GC12_perc ~ Normalised_Fourfold + Size)
          summary(mod_GC12_4x)
          
          # GC3 X Normalised_Fourfold
          ggplot(GC_gene, aes(x=GC3_perc,y=Normalised_Fourfold)) + geom_point() + geom_smooth(method = "lm") + theme_minimal() +
            ylab("Normalised Fourfold mutations") + xlab("GC3")
          # Test correlation
          cor.test(GC_gene$GC3_perc,GC_gene$Normalised_Fourfold)
          mod_GC3_4x <- lm(data = GC_gene, GC3_perc ~ Normalised_Fourfold + Size)
          summary(mod_GC3_4x)
        }
      }

    }
    }
    
}








