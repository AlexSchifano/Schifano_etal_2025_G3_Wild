setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")

library(seqinr)
library(stringr)

##### Comparing ancestral (Anc) and recent (Var) variants #####
{
  MEGA_CDS <- read.table("R_txt_files/Mutation_CDS.txt", header = T, sep="\t")
  MEGA_CDS_Syn <- MEGA_CDS[MEGA_CDS$Synonymous == "Synonymous",]
  Var_old <- read.table("R_txt_files/Variants_List.txt", header = T,sep="\t")

  ## Need to add a "Codon_From" Column for easier processing later in the loop
  for (ro in 1:nrow(MEGA_CDS_Syn)){
    # Extract to nucleotide from
    nuc_from <- substr(MEGA_CDS_Syn[ro,"Mutation"],1,1)
    # Extract the position in the codon
    cod_pos <- MEGA_CDS_Syn[ro,"Codon_Pos"]
    # Extract the Codon
    cod_to <- MEGA_CDS_Syn[ro,"Codon"]
    # Replace to recreate the Codon_from
    cod_from <- cod_to
    substr(cod_from,cod_pos,cod_pos) <- nuc_from
    # Add to the MEGA table
    MEGA_CDS_Syn[ro,"Codon_from"] <- cod_from
  }
  
  # Re-filter MEGA to have the true table
  Anc <- MEGA_CDS_Syn[MEGA_CDS_Syn$Strain_Name %in% Var_old$Strain_Name == F,]
  Var <- MEGA_CDS_Syn[MEGA_CDS_Syn$Strain_Name %in% Var_old$Strain_Name,]
  
  # Checking what is the prefered codon for each AA
  Codon_count_WS283 <- read.table("R_txt_files/Codon_counts_CDS_WS283.txt",header = T,sep="\t")
  # Separate the 2 serine codons
  Codon_count_WS283$AA[Codon_count_WS283$Codon %in% c("agt","agc","aga","agg")] <- "S1"
  Codon_count_WS283$AA[Codon_count_WS283$Codon %in% c("tct","tcc","tca","tcg")] <- "S2"
  AA_Anc_Var_Pref_count <- data.frame(AA = unique(Codon_count_WS283$AA),Prefered = "", To_Pref_Anc = 0, From_Pref_Anc = 0, To_Pref_Var = 0, From_Pref_Var = 0)
  for (AA in unique(Codon_count_WS283$AA)){
    
    AA_Count <- Codon_count_WS283[Codon_count_WS283$AA == AA,]
    Pref_codon <- AA_Count$Codon[AA_Count$Count == max(AA_Count$Count)]
    Pref_codon <- toupper(Pref_codon)
    AA_Anc_Var_Pref_count$Prefered[AA_Anc_Var_Pref_count$AA == AA] <- Pref_codon
    
    ## Check which muts are To/From the preferred codon
    AA_Anc_Var_Pref_count$To_Pref_Anc[AA_Anc_Var_Pref_count$AA == AA] <- nrow(Anc[Anc$Codon == Pref_codon,])
    AA_Anc_Var_Pref_count$To_Pref_Var[AA_Anc_Var_Pref_count$AA == AA] <- nrow(Var[Var$Codon == Pref_codon,])
    AA_Anc_Var_Pref_count$From_Pref_Anc[AA_Anc_Var_Pref_count$AA == AA] <- nrow(Anc[Anc$Codon_from == Pref_codon,])
    AA_Anc_Var_Pref_count$From_Pref_Var[AA_Anc_Var_Pref_count$AA == AA] <- nrow(Var[Var$Codon_from == Pref_codon,])
  }
  
  ## Stats
  {
    library(AMR)
    for(ro in 1:nrow(AA_Anc_Var_Pref_count)){
      
      a <- AA_Anc_Var_Pref_count[ro,"To_Pref_Anc"]
      b <- AA_Anc_Var_Pref_count[ro,"From_Pref_Anc"]
      c <- AA_Anc_Var_Pref_count[ro,"To_Pref_Var"]
      d <- AA_Anc_Var_Pref_count[ro,"From_Pref_Var"]
      
      mat <- matrix(c(a, b, c, d), nrow=2, byrow=TRUE)
      G <- g.test(mat)
      Fisher <- fisher.test(mat)
      
      AA_Anc_Var_Pref_count[ro,"G-Stat"] <- G$statistic
      AA_Anc_Var_Pref_count[ro,"G-pVal"] <- G$p.value
      AA_Anc_Var_Pref_count[ro,"Fisher_pVal"] <- Fisher$p.value
    }
  }
  
 }

##### Counting codons in the REF #####
{
  # Read the protein coding genome of the REF
  WS283_CDS <- read.fasta("../01-DATA/0-REF/FULL_REF/WS283_MtDNA_CDS_String.fasta")
  WS283_CDS_string <- paste(WS283_CDS[[1]],collapse = "")
  # splitting by codon
  WS283_CDS_Codon_split <- gsub("(.{3})", "\\1 ", WS283_CDS_string)
  WS283_CDS_Codon_split <- str_split(WS283_CDS_Codon_split,pattern = " ")
  
  ## Codon counts
  Codon_count_WS283 <- data.frame(Codon = unique(WS283_CDS_Codon_split[[1]]))
  # Remove empty
  Codon_count_WS283 <- data.frame(Codon = Codon_count_WS283[nchar(Codon_count_WS283$Codon) == 3,])
  # Count the number of each codon
  for (ro in 1:nrow(Codon_count_WS283)){
    codon <- Codon_count_WS283[ro,"Codon"]
    Codon_count_WS283[ro,"Count"] <- length(which(WS283_CDS_Codon_split[[1]] == codon))
    # Translate the codon
    codon <- as.vector(str_split(Codon_count_WS283[ro,"Codon"],pattern = ""))
    Codon_count_WS283[ro,"AA"] <- seqinr::translate(codon[[1]],numcode = 5)
  }

  # Write that table
  write.table(Codon_count_WS283, file = "R_txt_files/Codon_counts_CDS_WS283.txt", sep="\t", quote = F, row.names = F)
}

##### Relative observed freq of each codon #####
{
  Codon_count_WS283 <- read.table("R_txt_files/Codon_counts_CDS_WS283.txt",header = T,sep="\t")
  
  for (ro in 1:nrow(Codon_count_WS283)) {
    # Get the codon
    codonn <- Codon_count_WS283[ro,"Codon"]
    # The corresponding AA
    A <- Codon_count_WS283[ro,"AA"]
    # Add the counts for this AA
    tot_A <- sum(Codon_count_WS283$Count[Codon_count_WS283$AA == A])
    # Compute RSCU
    Codon_count_WS283[ro,"Obs_Freq"] <- Codon_count_WS283[ro,"Count"] / tot_A
  }
  
  # Write file
  write.table(Codon_count_WS283, file = "R_txt_files/Codon_counts_CDS_WS283.txt", sep="\t", quote = F, row.names = F)
}

##### Stat-tests to compare base composition at codon 3rd base #####
{
  # Get the table
  Codon_count_WS283 <- read.table("R_txt_files/Codon_counts_CDS_WS283.txt",header = T,sep="\t")
  
  ## Include still S but separating the two codon types
  Codon_count_WS283[Codon_count_WS283$AA == "S",]
  # I can see that rows 9,10,12 and 49 are one codon
  # Rows 19, 25, 47, 63 another type
  Codon_count_WS283[c(9,10,12,49),"AA"] <- "S1"
  Codon_count_WS283[c(19,25,47,63),"AA"] <- "S2"
  # I have to re-compute the observed freq to separate those
  Codon_count_WS283$Obs_Freq[Codon_count_WS283$AA == "S1"] <- Codon_count_WS283$Count[Codon_count_WS283$AA == "S1"] / 
    sum(Codon_count_WS283$Count[Codon_count_WS283$AA == "S1"])
  
  Codon_count_WS283$Obs_Freq[Codon_count_WS283$AA == "S2"] <- Codon_count_WS283$Count[Codon_count_WS283$AA == "S2"] / 
    sum(Codon_count_WS283$Count[Codon_count_WS283$AA == "S2"])
  
  AA_list <- data.frame(table(Codon_count_WS283$AA))
  
  #### At 4x codons ####
  {
    # Subset the table
    AA_list_4x <- AA_list[AA_list$Freq == 4,]
    Codon_4x <- Codon_count_WS283[Codon_count_WS283$AA %in% AA_list_4x$Var1,]
    
    # Get the third base
    Codon_4x$Base_3 <- substring(Codon_4x$Codon,3)
    
    # Make a table that counts for each 4x AA the number of base 3
    AA_Base_4x <- data.frame(AA=unique(Codon_4x$AA))
    row.names(AA_Base_4x) <- AA_Base_4x$AA
    for (ro in 1:nrow(AA_Base_4x)){
      amin <- AA_Base_4x[ro,"AA"]
      # sub table
      subtab <- Codon_4x[Codon_4x$AA == amin,]
      AA_Base_4x[ro,"T_count"] <- subtab$Count[subtab$Base_3 == "t"]
      AA_Base_4x[ro,"A_count"] <- subtab$Count[subtab$Base_3 == "a"]
      AA_Base_4x[ro,"G_count"] <- subtab$Count[subtab$Base_3 == "g"]
      AA_Base_4x[ro,"C_count"] <- subtab$Count[subtab$Base_3 == "c"]
    }
    
    # Compute total
    Total <- colSums(AA_Base_4x[,2:5])
    Total_freq <- Total / sum(Total)
    
    ## Overall g-test
    # Compare each AA
    # I remove row 2 because count < 5
    obs <- AA_Base_4x[,-1]
    expec <- Total_freq
    chisq.test(obs, p=expec_prop, simulate.p.value = T, B = 10000)
    
    # # Now I do ad-hoc tests to compare each amino acid counts to the counts of all other amino acids
    for (ro in 1:nrow(AA_Base_4x)) {
      amin <- AA_Base_4x[ro,"AA"]
      counts <- obs[amin,]
      colrm <- which(row.names(obs) == amin)
      expec <- colSums(obs[-colrm,])
      expec_prop <- expec / sum(expec)
      print(amin)
      print(chisq.test(counts, p = expec_prop, simulate.p.value = T, B = 10000))
    }
  }
  
  #### At 2x codons ####
  {
    # Subset the table
    AA_list_2x <- AA_list[AA_list$Freq == 2,]
    Codon_2x <- Codon_count_WS283[Codon_count_WS283$AA %in% AA_list_2x$Var1,]
    
    # Get the third base
    Codon_2x$Base_3 <- substring(Codon_2x$Codon,3)
    
    # Split the T|C ending and A|G ending codons
    Codon_2x_TC <- rbind(Codon_2x[Codon_2x$Base_3 == "c",],Codon_2x[Codon_2x$Base_3 == "t",])
    Codon_2x_AG <- rbind(Codon_2x[Codon_2x$Base_3 == "a",],Codon_2x[Codon_2x$Base_3 == "g",])
    
    ## Make a table that counts for each 2x AA in TC the number of base 3
    AA_Base_2x_TC <- data.frame(AA=unique(Codon_2x_TC$AA))
    row.names(AA_Base_2x_TC) <- AA_Base_2x_TC$AA
    for (ro in 1:nrow(AA_Base_2x_TC)){
      amin <- AA_Base_2x_TC[ro,"AA"]
      # sub table
      subtab <- Codon_2x_TC[Codon_2x_TC$AA == amin,]
      AA_Base_2x_TC[ro,"T_count"] <- subtab$Count[subtab$Base_3 == "t"]
      AA_Base_2x_TC[ro,"C_count"] <- subtab$Count[subtab$Base_3 == "c"]
    }
    ## Make a table that counts for each 2x AA in AG the number of base 3
    AA_Base_2x_AG <- data.frame(AA=unique(Codon_2x_AG$AA))
    row.names(AA_Base_2x_AG) <- AA_Base_2x_AG$AA
    for (ro in 1:nrow(AA_Base_2x_AG)){
      amin <- AA_Base_2x_AG[ro,"AA"]
      # sub table
      subtab <- Codon_2x_AG[Codon_2x_AG$AA == amin,]
      AA_Base_2x_AG[ro,"A_count"] <- subtab$Count[subtab$Base_3 == "a"]
      AA_Base_2x_AG[ro,"G_count"] <- subtab$Count[subtab$Base_3 == "g"]
    }
    
    # Compute total
    Total_TC <- colSums(AA_Base_2x_TC[,2:3])
    Total_freq_TC <- Total_TC / sum(Total_TC)
    
    Total_AG <- colSums(AA_Base_2x_AG[,2:3])
    Total_freq_AG <- Total_AG / sum(Total_AG)
    
    
    ## Overall g-test
    # Compare each AA
    obs_TC <- AA_Base_2x_TC[,-1]
    expec_TC <- Total_freq_TC
    chisq.test(obs_TC, p=expec_TC, simulate.p.value = T, B = 10000)
    ## Overall g-test
    # Compare each AA
    obs_AG <- AA_Base_2x_AG[,-1]
    expec_AG <- Total_freq_AG
    chisq.test(obs_AG, p=expec_AG, simulate.p.value = T, B = 10000)
    
    # # Now I do ad-hoc tests to compare each amino acid counts to the counts of all other amino acids
    for (ro in 1:nrow(AA_Base_2x_TC)) {
      amin <- AA_Base_2x_TC[ro,"AA"]
      counts <- obs_TC[amin,]
      colrm <- which(row.names(obs_TC) == amin)
      expec <- colSums(obs_TC[-colrm,])
      expec_prop <- expec / sum(expec)
      print(amin)
      print(chisq.test(counts, p = expec_prop, simulate.p.value = T, B = 10000))
    }
    # # Now I do ad-hoc tests to compare each amino acid counts to the counts of all other amino acids
    for (ro in 1:nrow(AA_Base_2x_AG)) {
      amin <- AA_Base_2x_AG[ro,"AA"]
      counts <- obs_AG[amin,]
      colrm <- which(row.names(obs_AG) == amin)
      expec <- colSums(obs_TC[-colrm,])
      expec_prop <- expec / sum(expec)
      print(amin)
      print(chisq.test(counts, p = expec_prop, simulate.p.value = T, B = 10000))
    }
  }
  
}

##### Comparing base composition at 2x and 4x #### 
{
  # Making the tables
  AA_list <- data.frame(table(Codon_count_WS283$AA))
  AA_list_4x <- AA_list[AA_list$Freq == 4,]
  
  Codon_count_WS283_4x <- Codon_count_WS283[Codon_count_WS283$AA %in% AA_list_4x$Var1,]
  Codon_count_WS283_2x <- Codon_count_WS283[Codon_count_WS283$AA %in% AA_list_2x$Var1,]
  
  # Add a column with the third base
  Codon_count_WS283_4x$Base_3 <- substring(Codon_count_WS283_4x$Codon,3,3)
  Codon_count_WS283_2x$Base_3 <- substring(Codon_count_WS283_2x$Codon,3,3)
  
  # Counting bases at 4x
  Base_comp_4x <- data.frame(A = sum(Codon_count_WS283_4x$Count[Codon_count_WS283_4x$Base_3 == "a"]),
                             T = sum(Codon_count_WS283_4x$Count[Codon_count_WS283_4x$Base_3 == "t"]),
                             C = sum(Codon_count_WS283_4x$Count[Codon_count_WS283_4x$Base_3 == "c"]),
                             G = sum(Codon_count_WS283_4x$Count[Codon_count_WS283_4x$Base_3 == "g"]))
  # Turn into %
  Base_comp_4x[2,] <- Base_comp_4x[1,] / sum(Base_comp_4x[1,])
  
  # Counting bases at 2x
  Base_comp_2x <- data.frame(A = sum(Codon_count_WS283_2x$Count[Codon_count_WS283_2x$Base_3 == "a"]),
                             T = sum(Codon_count_WS283_2x$Count[Codon_count_WS283_2x$Base_3 == "t"]),
                             C = sum(Codon_count_WS283_2x$Count[Codon_count_WS283_2x$Base_3 == "c"]),
                             G = sum(Codon_count_WS283_2x$Count[Codon_count_WS283_2x$Base_3 == "g"]))
  # Turn into %
  Base_comp_2x[2,] <- Base_comp_2x[1,] / sum(Base_comp_2x[1,])
  
  ## Stat test
  matrix <- matrix(as.numeric(c(Base_comp_4x[1,],Base_comp_2x[1,])),ncol = 4,byrow=T)
  library(AMR)
  g.test(matrix)
  chisq.test(matrix)
}

##### AT content simulation #####
{
  # i = AT content at 4x sites
  # j = GC content at 4x sites
  # i0 = initial AT content
  # j0 = initial GC content
  # a = rate of AT loss (A:T -> G|C:C|G)
  # b = rate of GC loss (C:G -> A|T:T|A)
  
  i0 = 667
  j0 = 112
  a = 0.56
  b = 0.44
  AT_evo <- data.frame(AT = i0, GC=j0, GEN = 0)
  
  # Loop for simulation
  for (gen in 2:100) {
    ix <- AT_evo[gen-1,1] - (AT_evo[gen-1,1] * a) + (AT_evo[gen-1,2] * b)
    jx <- AT_evo[gen-1,2] - (AT_evo[gen-1,2] * b) + (AT_evo[gen-1,1] * a)
    AT_evo[gen,] <- c(ix,jx,gen)
  }

  # Here we have that the equilibrium frequencies are ix / (ix+jx) = b and jx / (ix+jx) = a
  # At equilibrium, f(AC) = f(AT gain) and f(GC) = f(GC gain)
  eq_f(GC) = 0.56
  eq_f(AT) = 0.44
  
  obs_f(GC) = 0.14
  obs_f(AT) = 0.86
}


##### Expected freq of each codon under equal base use ##### WRONG
{
  # Get the table
  Codon_count_WS283 <- read.table("R_txt_files/Codon_counts_CDS_WS283.txt",header = T,sep="\t")
  
  ## Include still S but separating the two codon types
  Codon_count_WS283[Codon_count_WS283$AA == "S",]
  # I can see that rows 9,10,12 and 49 are one codon
  # Rows 19, 25, 47, 63 another type
  Codon_count_WS283[c(9,10,12,49),"AA"] <- "S1"
  Codon_count_WS283[c(19,25,47,63),"AA"] <- "S2"
  # I have to re-compute the observed freq to separate those
  Codon_count_WS283$Obs_Freq[Codon_count_WS283$AA == "S1"] <- Codon_count_WS283$Count[Codon_count_WS283$AA == "S1"] / 
    sum(Codon_count_WS283$Count[Codon_count_WS283$AA == "S1"])
  
  Codon_count_WS283$Obs_Freq[Codon_count_WS283$AA == "S2"] <- Codon_count_WS283$Count[Codon_count_WS283$AA == "S2"] / 
    sum(Codon_count_WS283$Count[Codon_count_WS283$AA == "S2"])
  
  AA_list <- data.frame(table(Codon_count_WS283$AA))
  
  #### Four-fold degenerate codons ####
  {
    # Subset the table to only the fourfold degenerate codons at 3rd base
    AA_list_4x <- AA_list[AA_list$Freq == 4,]
    Codon_count_WS283_4x <- Codon_count_WS283[Codon_count_WS283$AA %in% AA_list_4x$Var1,]
    
    # Add a column with the third base
    Codon_count_WS283_4x$Base_3 <- substring(Codon_count_WS283_4x$Codon,3,3)
    
    # Compute base composition at those codons' 3rd base
    Base_comp_4x <- data.frame(Nucleotide = c("a","c","g","t"), 
                               Count = c(sum(Codon_count_WS283_4x$Count[Codon_count_WS283_4x$Base_3 == "a"]),
                                         sum(Codon_count_WS283_4x$Count[Codon_count_WS283_4x$Base_3 == "c"]),
                                         sum(Codon_count_WS283_4x$Count[Codon_count_WS283_4x$Base_3 == "g"]),
                                         sum(Codon_count_WS283_4x$Count[Codon_count_WS283_4x$Base_3 == "t"])))
    Base_comp_4x$Freq <- Base_comp_4x$Count / sum(Base_comp_4x$Count)
    
    # Loop to compute expected counts and freq under equal base usage
    for (ro in 1:nrow(Codon_count_WS283_4x)) {
      codon <- Codon_count_WS283_4x[ro,"Codon"]
      AA <- Codon_count_WS283_4x[ro,"AA"]
      # Get the 3rd base
      Base_3 <- substring(codon,3)
      # Get the total AA count
      AA_count <- sum(Codon_count_WS283$Count[Codon_count_WS283$AA == AA])
      # Putting the expected freq (equal freq of base at 4x)
      Codon_count_WS283_4x[ro,"Expected_freq_equal_base"] <- Base_comp_4x$Freq[Base_comp_4x$Nucleotide == Base_3]
      # Expected Count
      Codon_count_WS283_4x[ro,"Expected_count_equal_base"] <- Codon_count_WS283_4x[ro,"Expected_freq_equal_base"]*AA_count
    }
    
    # NB_RSCU as the ratio observed/expected
    Codon_count_WS283_4x$NB_RSCU <- Codon_count_WS283_4x$Obs_Freq / Codon_count_WS283_4x$Expected_freq_equal_base
    
    # Plot the relationship between Codon count and NB_RSCU
    library(ggplot2)
    CUB_4x <- ggplot(Codon_count_WS283_4x,aes(x=Count,y=NB_RSCU,label=AA,col=Base_3)) + geom_point() + geom_label() + 
      ylim(c(0,2.5)) + scale_colour_manual(values = c("a" = "red","c" = "blue","g" = "turquoise3","t" = "orange")) +
      theme_minimal()
    
    ggsave(filename = "Writing_NEW/5-Figures/CUB/CUB_4x.pdf",width=20,height = 13,plot = CUB_4x)
    
    # Order the table 
    Codon_count_WS283_4x <- Codon_count_WS283_4x[order(Codon_count_WS283_4x$Codon),]
    
    # Writing the table
    write.table(Codon_count_WS283_4x, file = "R_txt_files/Codon_NBRSCU_Equal_Base_4x.txt", quote = F, row.names = F, sep = "\t")
    
    ## Fisher test on this data to compare observed and expected
    for (amin in unique(Codon_count_WS283_4x$AA)) {
      obs <- Codon_count_WS283_4x$Count[Codon_count_WS283_4x$AA == amin]
      exp <- Codon_count_WS283_4x$Expected_count_equal_base[Codon_count_WS283_4x$AA == amin]
      exp <- round(exp)
      matrix <- matrix(c(obs,exp),byrow = T,ncol = 4)
      print(amin)
      print(fisher.test(matrix)) 
    }
    
  }
  
  #### Two-fold degenerate sites ####
  {
    # Subset the table
    AA_list_2x <- AA_list[AA_list$Freq == 2,]
    Codon_count_WS283_2x <- Codon_count_WS283[Codon_count_WS283$AA %in% AA_list_2x$Var1,]
    
    # Base composition at those sites
    Codon_count_WS283_2x$Base_3 <- substring(Codon_count_WS283_2x$Codon,3,3)
    Base_comp_2x <- data.frame(a = sum(Codon_count_WS283_2x$Count[Codon_count_WS283_2x$Base_3 == "a"]),
                               t = sum(Codon_count_WS283_2x$Count[Codon_count_WS283_2x$Base_3 == "t"]),
                               c = sum(Codon_count_WS283_2x$Count[Codon_count_WS283_2x$Base_3 == "c"]),
                               g = sum(Codon_count_WS283_2x$Count[Codon_count_WS283_2x$Base_3 == "g"]))
    Base_comp_2x[2,] <- Base_comp_2x[1,] / sum(Base_comp_2x[1,])
    row.names(Base_comp_2x) <- c("Count","Proportion")
    
    for (ro in 1:nrow(Codon_count_WS283_2x)){
      amin <- Codon_count_WS283_2x[ro,"AA"]
      base <- Codon_count_WS283_2x[ro,"Base_3"]
      # Subset the table to only different to the base
      nobase <- Codon_count_WS283_2x[Codon_count_WS283_2x$Base_3 != base,]
      # Get the second base pair that leads to the same codon
      base_2 <- nobase$Base_3[nobase$AA == amin]
      # Compute the expected freq
      exp_freq_base <- Base_comp_2x["Count",base] / sum(Base_comp_2x["Count",base] + Base_comp_2x["Count",base_2])
      exp_freq_base_2 <- Base_comp_2x["Count",base_2] / sum(Base_comp_2x["Count",base] + Base_comp_2x["Count",base_2])
      # Put in the table
      Codon_count_WS283_2x[ro,"Expected_freq"] <- exp_freq_base
      # Put the expected counts
      count_amin <- sum(Codon_count_WS283_2x$Count[Codon_count_WS283_2x$AA == amin])
      Codon_count_WS283_2x[ro,"Expected_count"] <- exp_freq_base * count_amin
    }
    # Adding NB_RSCU
    Codon_count_WS283_2x$NB_RSCU <- Codon_count_WS283_2x$Count / Codon_count_WS283_2x$Expected_count
    
    # Plot the relationship between Codon count and NB_RSCU
    library(ggplot2)
    CUB_2x <- ggplot(Codon_count_WS283_2x,aes(x=Count,y=NB_RSCU,label=AA,col=Base_3)) + geom_point() + geom_label() + 
      ylim(c(0,2.5)) + scale_colour_manual(values = c("a" = "red","c" = "blue","g" = "turquoise3","t" = "orange")) +
      theme_minimal()
    
    ggsave(filename = "Writing_NEW/5-Figures/CUB/CUB_2x.pdf",width=20,height = 13,plot = CUB_2x)
    
    # Order the table 
    Codon_count_WS283_2x <- Codon_count_WS283_2x[order(Codon_count_WS283_2x$Codon),]
    Codon_count_WS283_2x <- Codon_count_WS283_2x[order(Codon_count_WS283_2x$AA),]
    
    # Writing the table
    write.table(Codon_count_WS283_2x, file = "R_txt_files/Codon_NBRSCU_Equal_Base_4x.txt", quote = F, row.names = F, sep = "\t")
    
    ## Fisher test on this data to compare observed and expected
    for (amin in unique(Codon_count_WS283_2x$AA)) {
      obs <- Codon_count_WS283_2x$Count[Codon_count_WS283_2x$AA == amin]
      exp <- Codon_count_WS283_2x$Expected_count[Codon_count_WS283_2x$AA == amin]
      exp <- round(exp)
      matrix <- matrix(c(obs,exp),byrow = T,ncol = 2)
      print(amin)
      print(fisher.test(matrix)) 
    }
  }
}
