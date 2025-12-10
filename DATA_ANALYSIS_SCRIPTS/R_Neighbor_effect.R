setwd("~/Documents/C.elegans_sync/Final_files/")

library(seqinr)
### Investigate the neighbor effect

############### Making the neighbor effect table ###############
{
  # I look at each mutation and need to record the upstream and downstream nucleotide.
  # Read mut table, only syn mutations
  MEGA <- read.table("02-ANALYSES/R_txt_files/Mutation_List.txt", header = T,sep="\t")
  MEGA <- MEGA[MEGA$Synonymous == "Synonymous",]
  
  # Remove start and stop codon mutations
  GTF <- read.table("01-DATA/0-REF/WS283_MtDNA.gtf")
  GTF_CDS <- GTF[GTF$V5 == "CDS",]
  MEGA <- MEGA[MEGA$POS_REF %in% GTF_CDS$V3 == F,]
  MEGA <- MEGA[MEGA$POS_REF %in% GTF_CDS$V2 == F,]
  
  # I need to split it because the fasta files are not the same depending on if internal or leaf mutation
  Anc_list <- read.table("01-DATA/8-COUNT_MUTATIONS/List_Anc.txt")
  MEGA_Str <- MEGA[(MEGA$Strain_Name %in% Anc_list$V1)== F,]
  MEGA_Anc <- MEGA[MEGA$Strain_Name %in% Anc_list$V1,]
  
  # Table for storage
  Neigh_Str <- data.frame(row.names = row.names(MEGA_Str))
  Neigh_Anc <- data.frame(row.names = row.names(MEGA_Anc))
  
  ## Loop through rows of MEGA_Str
  for (ro in 1:nrow(MEGA_Str)) {
    # store position
    pos <- MEGA_Str[ro,"POS_Full"]
    # strain
    strain <- MEGA_Str[ro,"Ancestor"]
    # retrieve the fasta
    fasta_name <- paste("01-DATA/9-VARIANT_CALLING/REF/",strain,".fa",sep="")
    fasta <- read.fasta(fasta_name,forceDNAtolower = F)
    # retrieve the upstream and downstream bp
    pos_nuc <- fasta[[1]][pos]
    pos_up <- fasta[[1]][pos-1]
    pos_dw <- fasta[[1]][pos+1]
    # add info to table
    Neigh_Str[ro,"POS_Full"] <- pos
    Neigh_Str[ro,"POS_REF"] <- MEGA_Str[ro,"POS_REF"]
    Neigh_Str[ro,"Mutation"] <- MEGA_Str[ro,"Mutation"]
    Neigh_Str[ro,"Strain"] <- MEGA_Str[ro,"Strain_Name"]
    Neigh_Str[ro,"UP"] <- pos_up
    Neigh_Str[ro,"DOWN"] <- pos_dw
    Neigh_Str[ro,"Degeneracy"] <- MEGA_Str[ro,"Site_4x_2x"]
    Neigh_Str[ro,"Codon_pos"] <- MEGA_Str[ro,"Codon_Pos"]
  }
  
  ## Loop through rows of MEGA_Anc
  # I need to extract the sequences from the MSA
  MSA <- read.fasta("01-DATA/8-COUNT_MUTATIONS/Anc.fasta",forceDNAtolower = F)
  for (ro in 1:nrow(MEGA_Anc)) {
    # store position
    pos <- MEGA_Anc[ro,"POS_Full"]
    # strain
    strain <- MEGA_Anc[ro,"Strain_Name"]
    # retrieve the fasta
    fasta <- MSA[strain]
    # retrieve the upstream and downstream bp
    pos_nuc <- fasta[[1]][pos]
    pos_up <- fasta[[1]][pos-1]
    pos_dw <- fasta[[1]][pos+1]
    # add info to table
    Neigh_Anc[ro,"POS_Full"] <- pos
    Neigh_Anc[ro,"POS_REF"] <- MEGA_Anc[ro,"POS_REF"]
    Neigh_Anc[ro,"Mutation"] <- MEGA_Anc[ro,"Mutation"]
    Neigh_Anc[ro,"Strain"] <- strain
    Neigh_Anc[ro,"UP"] <- pos_up
    Neigh_Anc[ro,"DOWN"] <- pos_dw
    Neigh_Anc[ro,"Degeneracy"] <- MEGA_Anc[ro,"Site_4x_2x"]
    Neigh_Anc[ro,"Codon_pos"] <- MEGA_Anc[ro,"Codon_Pos"]
  }
  
  ## Save the files
  # Combine first
  Neigh_Anc$Strain_type <- "Anc"
  Neigh_Str$Strain_type <- "Var"
  Neighbor <- rbind(Neigh_Anc,Neigh_Str)
  write.table(Neighbor, file = "02-ANALYSES/R_txt_files/Neighbor_list.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  library(openxlsx)
  write.xlsx(Neighbor, file = "02-ANALYSES/R_txt_files/Neighbor_list.xlsx")
  
}

###############  Now that I have the list I can check for which upstream / downstream nucleotide is more often linked with mutation ############### 
{
  Neighbor <- read.table("02-ANALYSES/R_txt_files/Neighbor_list.txt", header = T, sep = "\t")
  Neigh_Anc <- Neighbor[Neighbor$Strain_type == "Anc",]
  Neigh_Str <- Neighbor[Neighbor$Strain_type == "Var",]
  
  ############# First look ############# 
  {
    # Raw proportions 
    raw_prop_anc_up <- table(Neigh_Anc$UP) / sum(table(Neigh_Anc$UP))
    raw_prop_str_up <- table(Neigh_Str$UP) / sum(table(Neigh_Str$UP))
    
    raw_prop_anc_dw <- table(Neigh_Anc$DOWN) / sum(table(Neigh_Anc$DOWN))
    raw_prop_str_dw <- table(Neigh_Str$DOWN) / sum(table(Neigh_Str$DOWN))
    
    # Base composition
    ws283 <- read.fasta("01-DATA/0-REF/WS283_MtDNA.fasta",forceDNAtolower = F)
    bases <- seqinr::count(ws283[[1]], wordsize = 1, alphabet = c("A","T","C","G")) 
    bases_prop <- bases/sum(bases)
    
    # Normalised
    raw_prop_anc_up / bases_prop 
    raw_prop_str_up / bases_prop 
    raw_prop_anc_dw / bases_prop 
    raw_prop_str_dw / bases_prop 
  }
  
  ############# Looking at 4x 2x separately ############# 
  {
    ###### Load the script with the functions I made ######
    source("~/Documents/Scripts/R_FUNCTIONS.R")
    
    ##### Degeneracy tag of ref #####
    {
      ## Base composition of only CDS
      REF <- read.fasta("01-DATA/0-REF/REF_CDS_Conca.fasta",forceDNAtolower = F)
      degeneracy_tag_CDS(REF, Genetic_code = 5)

      ######## Compute base compositions
      # Upstream and downstream of 4x pos
      {
        ## The list of positions to screen (all position except for START and STOP)
        pos_screen <- c(1:nrow(Seq_table))
        # I need the start and stop pos of the genes in the CDS concatenated to remove them (not 2x 4x behavior) + only 1 mut
        start_stop_pos <- c(1,435,436,669,670,1545,1546,2145,2146,2994,2995,4107,4108,4875,4876,6105,6106,7683,7684,8379,8380,8715,8716,10299)
        pos_screen <- pos_screen[-start_stop_pos]
        
        Neighbor_base_composition_REF <- data.frame(row.names = c("A","T","C","G"), UP_2x = rep(0,4), DOWN_2x = 0, UP_4x = 0, DOWN_4x = 0)
        for (pos in pos_screen){
          if (Seq_table[pos,"Degeneracy"] == "Twofold"){
            pos_up <- pos - 1
            pos_dw <- pos + 1
            base_up <- Seq_table[pos_up,"Seq"]
            base_dw <- Seq_table[pos_dw,"Seq"]
            Neighbor_base_composition_REF[base_up,"UP_2x"] <- Neighbor_base_composition_REF[base_up,"UP_2x"] + 1
            Neighbor_base_composition_REF[base_dw,"DOWN_2x"] <- Neighbor_base_composition_REF[base_dw,"DOWN_2x"] + 1
          }
          if (Seq_table[pos,"Degeneracy"] == "Fourfold"){
            pos_up <- pos - 1
            pos_dw <- pos + 1
            base_up <- Seq_table[pos_up,"Seq"]
            base_dw <- Seq_table[pos_dw,"Seq"]
            Neighbor_base_composition_REF[base_up,"UP_4x"] <- Neighbor_base_composition_REF[base_up,"UP_4x"] + 1
            Neighbor_base_composition_REF[base_dw,"DOWN_4x"] <- Neighbor_base_composition_REF[base_dw,"DOWN_4x"] + 1
        }
        }
        
        # Proportions
        {
          Prop_Neighbor_base_composition_REF <- Neighbor_base_composition_REF
          Prop_Neighbor_base_composition_REF$UP_2x <- Prop_Neighbor_base_composition_REF$UP_2x / sum(Prop_Neighbor_base_composition_REF$UP_2x)
          Prop_Neighbor_base_composition_REF$UP_4x <- Prop_Neighbor_base_composition_REF$UP_4x / sum(Prop_Neighbor_base_composition_REF$UP_4x)
          Prop_Neighbor_base_composition_REF$DOWN_2x <- Prop_Neighbor_base_composition_REF$DOWN_2x / sum(Prop_Neighbor_base_composition_REF$DOWN_2x)
          Prop_Neighbor_base_composition_REF$DOWN_4x <- Prop_Neighbor_base_composition_REF$DOWN_4x / sum(Prop_Neighbor_base_composition_REF$DOWN_4x)
        }

      }
      }

    ### Fourfold ==> I can normalize by base composition up and down stream of 4x mutations
    {
      Neighbor_4x <- Neighbor[Neighbor$Degeneracy == "4x",]
      Neighbor_count_4x <- data.frame(row.names = row.names(Neighbor_base_composition_REF), UP = rep(0,4), DOWN = rep(0,4))
      for (bp in row.names(Neighbor_count_4x)) {
        print(bp)
        Neighbor_count_4x[bp,"UP"] <- nrow(Neighbor_4x[Neighbor_4x$UP == bp,])
        Neighbor_count_4x[bp,"DOWN"] <- nrow(Neighbor_4x[Neighbor_4x$DOWN == bp,])
      }
      Neighbor_count_4x$Normalized_UP <- Neighbor_count_4x$UP / Neighbor_base_composition_REF$UP_4x
      Neighbor_count_4x$Normalized_DOWN <- Neighbor_count_4x$DOWN / Neighbor_base_composition_REF$DOWN_4x
    }

    ### Twofold ==> I can normalize by base composition up and down stream of 2x mutations
    {
      Neighbor_2x <- Neighbor[Neighbor$Degeneracy == "2x",]
      Neighbor_count_2x <- data.frame(row.names = row.names(Neighbor_base_composition_REF), UP = rep(0,4), DOWN = rep(0,4))
      for (bp in row.names(Neighbor_count_2x)) {
        print(bp)
        Neighbor_count_2x[bp,"UP"] <- nrow(Neighbor_2x[Neighbor_2x$UP == bp,])
        Neighbor_count_2x[bp,"DOWN"] <- nrow(Neighbor_2x[Neighbor_2x$DOWN == bp,])
      }
      Neighbor_count_2x$Normalized_UP <- Neighbor_count_2x$UP / Neighbor_base_composition_REF$UP_2x
      Neighbor_count_2x$Normalized_DOWN <- Neighbor_count_2x$DOWN / Neighbor_base_composition_REF$DOWN_2x
    }

    ## Making a final table
    {
      Neighbor_Count <- data.frame(row.names = row.names(Neighbor_base_composition_REF), Composition_2x_UP = Neighbor_base_composition_REF$UP_2x,
                                   Composition_4x_UP = Neighbor_base_composition_REF$UP_4x, Composition_2x_DOWN = Neighbor_base_composition_REF$DOWN_2x,
                                   Composition_4x_DOWN = Neighbor_base_composition_REF$DOWN_4x, Count_Mut_2x_UP = Neighbor_count_2x$UP, 
                                   Count_Mut_4x_UP = Neighbor_count_4x$UP, Count_Mut_2x_DOWN = Neighbor_count_2x$DOWN, Count_Mut_4x_DOWN = Neighbor_count_4x$DOWN,
                                   Normalised_Count_2x_UP = Neighbor_count_2x$Normalized_UP, Normalised_Count_4x_UP = Neighbor_count_4x$Normalized_UP, 
                                   Normalised_Count_2x_DOWN = Neighbor_count_2x$Normalized_DOWN, Normalised_Count_4x_DOWN = Neighbor_count_4x$Normalized_DOWN)
      # Divided by 0 fix
      Neighbor_Count["A","Normalised_Count_4x_UP"] <- 0
    }

    ## Plotting the distribution of neighbors
    {
      # UP 2x
      library(ggplot2)
      ggplot(Neighbor_Count,aes(x = row.names(Neighbor_Count), y = Normalised_Count_2x_UP)) + geom_point() + ggtitle(label = "Upstream of Twofold") +
        ylim(c(0,1)) + theme_minimal() + xlab("Nucleotide") + ylab("Normalised count (Count / Base composition)") 
      # UP 4x
      ggplot(Neighbor_Count,aes(x = row.names(Neighbor_Count), y = Normalised_Count_4x_UP)) + geom_point() + ggtitle(label = "Upstream of Fourfold") +
        ylim(c(0,1)) + theme_minimal() + xlab("Nucleotide") + ylab("Normalised count (Count / Base composition)") 
      # DOWN 2x
      ggplot(Neighbor_Count,aes(x = row.names(Neighbor_Count), y = Normalised_Count_2x_DOWN)) + geom_point() + ggtitle(label = "Downstream of Twofold") +
        ylim(c(0,2)) + theme_minimal() + xlab("Nucleotide") + ylab("Normalised count (Count / Base composition)") 
      # DOWN 4x
      ggplot(Neighbor_Count,aes(x = row.names(Neighbor_Count), y = Normalised_Count_4x_DOWN)) + geom_point() + ggtitle(label = "Downstream of Fourfold") +
        ylim(c(0,2)) + theme_minimal() + xlab("Nucleotide") + ylab("Normalised count (Count / Base composition)") 
      
      # Reshape data for ggplot
      Neighbor_Count_long <- data.frame(
        Nucleotide = rep(row.names(Neighbor_Count), times = 4),
        Count = c(Neighbor_Count$Normalised_Count_2x_UP, 
                  Neighbor_Count$Normalised_Count_2x_DOWN,
                  Neighbor_Count$Normalised_Count_4x_DOWN, 
                  Neighbor_Count$Normalised_Count_4x_UP),
        Direction = rep(c("Up", "Down", "Down", "Up"), each = nrow(Neighbor_Count)),
        Degeneracy = rep(c("2x", "2x", "4x", "4x"), each = nrow(Neighbor_Count))) # Changed from "Fold_Change"
      
      ggplot(Neighbor_Count_long, aes(x = Nucleotide, y = Count, shape = Direction, fill = Degeneracy)) +
        geom_point(size = 3, color = "black") +  # Black outline for visibility
        scale_shape_manual(values = c("Up" = 24, "Down" = 25)) +  # Shapes for Up/Down
        scale_fill_manual(values = c("2x" = "orange", "4x" = "turquoise3"), guide = guide_legend(override.aes = list(shape = 21, color = "black"))) +  
        theme_minimal() +
        xlab("Nucleotide") +
        ylab("Normalised count (Count / Base composition)") +
        labs(shape = "Direction", fill = "Degeneracy")  # Correct legend labels
    }
    
    ## Grouping AT GC
    {
      Grouped_Neighbor <- data.frame(Nucleotide = c("A|T","C|G","A|T","C|G","A|T","C|G","A|T","C|G"), 
                                     Degeneracy = c("Twofold","Twofold","Twofold","Twofold","Fourfold","Fourfold","Fourfold","Fourfold"),
                                     Direction = c("5'","5'","3'","3'","5'","5'","3'","3'"),
                                     Count = 0, Composition = 0,
                                     row.names = c("2x_UP_AT","2x_UP_GC","2x_DOWN_AT","2x_DOWN_GC","4x_UP_AT","4x_UP_GC","4x_DOWN_AT","4x_DOWN_GC"))
      ## Summing
      {
        # 2x UP
        Grouped_Neighbor["2x_UP_AT","Count"] <- Neighbor_Count["A","Count_Mut_2x_UP"] + Neighbor_Count["T","Count_Mut_2x_UP"]
        Grouped_Neighbor["2x_UP_AT","Composition"] <- Neighbor_Count["A","Composition_2x_UP"] + Neighbor_Count["T","Composition_2x_UP"]
        Grouped_Neighbor["2x_UP_GC","Count"] <- Neighbor_Count["G","Count_Mut_2x_UP"] + Neighbor_Count["C","Count_Mut_2x_UP"]
        Grouped_Neighbor["2x_UP_GC","Composition"] <- Neighbor_Count["G","Composition_2x_UP"] + Neighbor_Count["C","Composition_2x_UP"]
        # 2x DOWN
        Grouped_Neighbor["2x_DOWN_AT","Count"] <- Neighbor_Count["A","Count_Mut_2x_DOWN"] + Neighbor_Count["T","Count_Mut_2x_DOWN"]
        Grouped_Neighbor["2x_DOWN_AT","Composition"] <- Neighbor_Count["A","Composition_2x_DOWN"] + Neighbor_Count["T","Composition_2x_DOWN"]
        Grouped_Neighbor["2x_DOWN_GC","Count"] <- Neighbor_Count["G","Count_Mut_2x_DOWN"] + Neighbor_Count["C","Count_Mut_2x_DOWN"]
        Grouped_Neighbor["2x_DOWN_GC","Composition"] <- Neighbor_Count["G","Composition_2x_DOWN"] + Neighbor_Count["C","Composition_2x_DOWN"]
        # 4x UP
        Grouped_Neighbor["4x_UP_AT","Count"] <- Neighbor_Count["A","Count_Mut_4x_UP"] + Neighbor_Count["T","Count_Mut_4x_UP"]
        Grouped_Neighbor["4x_UP_AT","Composition"] <- Neighbor_Count["A","Composition_4x_UP"] + Neighbor_Count["T","Composition_4x_UP"]
        Grouped_Neighbor["4x_UP_GC","Count"] <- Neighbor_Count["G","Count_Mut_4x_UP"] + Neighbor_Count["C","Count_Mut_4x_UP"]
        Grouped_Neighbor["4x_UP_GC","Composition"] <- Neighbor_Count["G","Composition_4x_UP"] + Neighbor_Count["C","Composition_4x_UP"]
        # 4x DOWN
        Grouped_Neighbor["4x_DOWN_AT","Count"] <- Neighbor_Count["A","Count_Mut_4x_DOWN"] + Neighbor_Count["T","Count_Mut_4x_DOWN"]
        Grouped_Neighbor["4x_DOWN_AT","Composition"] <- Neighbor_Count["A","Composition_4x_DOWN"] + Neighbor_Count["T","Composition_4x_DOWN"]
        Grouped_Neighbor["4x_DOWN_GC","Count"] <- Neighbor_Count["G","Count_Mut_4x_DOWN"] + Neighbor_Count["C","Count_Mut_4x_DOWN"]
        Grouped_Neighbor["4x_DOWN_GC","Composition"] <- Neighbor_Count["G","Composition_4x_DOWN"] + Neighbor_Count["C","Composition_4x_DOWN"]
      }
      
      ## Normalized 
      Grouped_Neighbor$Normalized <- Grouped_Neighbor$Count / Grouped_Neighbor$Composition
      
      ## Plotting
      {
        ggplot(Grouped_Neighbor, aes(x = Nucleotide, y = Normalized, shape = Direction, fill = Degeneracy)) +
          geom_point(size = 3, color = "black") +  # Black outline for visibility
          scale_shape_manual(values = c("5'" = 24, "3'" = 25)) +  # Shapes for Up/Down
          scale_fill_manual(values = c("Twofold" = "orange", "Fourfold" = "turquoise3"), guide = guide_legend(override.aes = list(shape = 21, color = "black"))) +  
          theme_minimal() +
          xlab("Nucleotide") +
          ylab("Normalised count (Count / Base composition)") +
          labs(shape = "Direction", fill = "Degeneracy")  # Correct legend labels
        
        
        Neigh <- ggplot(Grouped_Neighbor, aes(x = Nucleotide, y = Normalized, shape = Direction, fill = Degeneracy)) +
          geom_point(size = 5, color = "black") +  # Black outline for visibility
          geom_line(aes(group = interaction(Direction, Degeneracy)), color = "red") +  # Add lines between points with the same Direction and Degeneracy
          scale_shape_manual(values = c("5'" = 24, "3'" = 25)) +  # Shapes for Up/Down
          scale_fill_manual(values = c("Twofold" = "orange", "Fourfold" = "turquoise3"), guide = guide_legend(override.aes = list(shape = 21, color = "black"))) +  
          theme_minimal() +
          xlab("Nucleotide") +
          ylab("Composition-normalized neighbour counts") +
          labs(shape = "Direction", fill = "Degeneracy")  # Correct legend labels
        
        ggsave(Neigh, filename = "03-SV_ANALYSES/Neighbour_Effect/Fig_Neighbour_effect.pdf")
      }

      ## Compute expected
      # Rate mut (DOWN and UP should be equal)
      {
        rate_2x_UP <- (sum(Grouped_Neighbor["2x_UP_AT","Count"], Grouped_Neighbor["2x_UP_GC","Count"]) /
                         sum(Grouped_Neighbor["2x_UP_AT","Composition"], Grouped_Neighbor["2x_UP_GC","Composition"]) )
        
        rate_2x_DOWN <- (sum(Grouped_Neighbor["2x_DOWN_AT","Count"], Grouped_Neighbor["2x_DOWN_GC","Count"]) /
                           sum(Grouped_Neighbor["2x_DOWN_AT","Composition"], Grouped_Neighbor["2x_DOWN_GC","Composition"]) )
        
        rate_4x_UP <- (sum(Grouped_Neighbor["4x_UP_AT","Count"], Grouped_Neighbor["4x_UP_GC","Count"]) /
                         sum(Grouped_Neighbor["4x_UP_AT","Composition"], Grouped_Neighbor["4x_UP_GC","Composition"]) )
        
        rate_4x_DOWN <- (sum(Grouped_Neighbor["4x_DOWN_AT","Count"], Grouped_Neighbor["4x_DOWN_GC","Count"]) /
                           sum(Grouped_Neighbor["4x_DOWN_AT","Composition"], Grouped_Neighbor["4x_DOWN_GC","Composition"]) )
      }
       
      # Expected count = rate * composition
      {
        Grouped_Neighbor["2x_UP_AT","Expected_count"] <- rate_2x_UP * Grouped_Neighbor["2x_UP_AT","Composition"]
        Grouped_Neighbor["2x_DOWN_AT","Expected_count"] <- rate_2x_DOWN * Grouped_Neighbor["2x_DOWN_AT","Composition"]
        Grouped_Neighbor["4x_UP_AT","Expected_count"] <- rate_4x_UP * Grouped_Neighbor["4x_UP_AT","Composition"]
        Grouped_Neighbor["4x_DOWN_AT","Expected_count"] <- rate_4x_DOWN * Grouped_Neighbor["4x_DOWN_AT","Composition"]
        Grouped_Neighbor["2x_UP_GC","Expected_count"] <- rate_2x_UP * Grouped_Neighbor["2x_UP_GC","Composition"]
        Grouped_Neighbor["2x_DOWN_GC","Expected_count"] <- rate_2x_DOWN * Grouped_Neighbor["2x_DOWN_GC","Composition"]
        Grouped_Neighbor["4x_UP_GC","Expected_count"] <- rate_4x_UP * Grouped_Neighbor["4x_UP_GC","Composition"]
        Grouped_Neighbor["4x_DOWN_GC","Expected_count"] <- rate_4x_DOWN * Grouped_Neighbor["4x_DOWN_GC","Composition"]
      }
      
      # Expected normalized count
      Grouped_Neighbor$Expected_Normalized_count <- Grouped_Neighbor$Expected_count / Grouped_Neighbor$Composition
      
      # G-test
      {
        # 2x Up
        g.test(c(Grouped_Neighbor["2x_UP_AT","Count"], Grouped_Neighbor["2x_UP_GC","Count"]), 
               p = c(Grouped_Neighbor["2x_UP_AT","Expected_count"]/sum(Grouped_Neighbor["2x_UP_AT","Expected_count"], Grouped_Neighbor["2x_UP_GC","Expected_count"]), 
                     Grouped_Neighbor["2x_UP_GC","Expected_count"]/sum(Grouped_Neighbor["2x_UP_AT","Expected_count"], Grouped_Neighbor["2x_UP_GC","Expected_count"])))
        # 2x Down
        g.test(c(Grouped_Neighbor["2x_DOWN_AT","Count"], Grouped_Neighbor["2x_DOWN_GC","Count"]), 
               p = c(Grouped_Neighbor["2x_DOWN_AT","Expected_count"]/sum(Grouped_Neighbor["2x_DOWN_AT","Expected_count"], Grouped_Neighbor["2x_DOWN_GC","Expected_count"]), 
                     Grouped_Neighbor["2x_DOWN_GC","Expected_count"]/sum(Grouped_Neighbor["2x_DOWN_AT","Expected_count"], Grouped_Neighbor["2x_DOWN_GC","Expected_count"])))
        
        # 4x Up
        g.test(c(Grouped_Neighbor["4x_UP_AT","Count"], Grouped_Neighbor["4x_UP_GC","Count"]), 
               p = c(Grouped_Neighbor["4x_UP_AT","Expected_count"]/sum(Grouped_Neighbor["4x_UP_AT","Expected_count"], Grouped_Neighbor["4x_UP_GC","Expected_count"]), 
                     Grouped_Neighbor["4x_UP_GC","Expected_count"]/sum(Grouped_Neighbor["4x_UP_AT","Expected_count"], Grouped_Neighbor["4x_UP_GC","Expected_count"])))
        
        # 4x Down
        g.test(c(Grouped_Neighbor["4x_DOWN_AT","Count"], Grouped_Neighbor["4x_DOWN_GC","Count"]), 
               p = c(Grouped_Neighbor["4x_DOWN_AT","Expected_count"]/sum(Grouped_Neighbor["4x_DOWN_AT","Expected_count"], Grouped_Neighbor["4x_DOWN_GC","Expected_count"]), 
                     Grouped_Neighbor["4x_DOWN_GC","Expected_count"]/sum(Grouped_Neighbor["4x_DOWN_AT","Expected_count"], Grouped_Neighbor["4x_DOWN_GC","Expected_count"])))
        
      }
      
      }
    
    
  }
}




