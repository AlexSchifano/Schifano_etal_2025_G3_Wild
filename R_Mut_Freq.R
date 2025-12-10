setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")

###### Load the script with the functions I made ######
source("~/Documents/Scripts/R_FUNCTIONS.R")

## Load mutation table
MEGA <- read.table("R_txt_files/Mutation_List.txt",header=T,sep="\t")

## Load GTF and name the columns
V170_REF_GTF <- read.table("../01-DATA/9-VARIANT_CALLING/GTF/GTF_V170.fa.gtf")
colnames(V170_REF_GTF) <- c("Chr","START","END","Gene_ID","Gene_Type")
# Need to add a column with size
V170_REF_GTF$SIZE <- V170_REF_GTF$END - V170_REF_GTF$START + 1

##### For all genes but not intergenic #####
mut_freq_gene(MEGA)

##### For CDS #####
{
  # Only getting CDS muts
  MEGA_CDS <- MEGA[MEGA$Gene_Type == "CDS",]
  # Only getting CDS from the GTF
  V170_REF_GTF_CDS <- V170_REF_GTF[V170_REF_GTF$Gene_Type == "CDS",]
  # Using the function to get gene mut counts and freq etc
  mut_freq_gene(MEGA_CDS,V170_REF_GTF_CDS)
}

##### For CDS Syn #####
{
  # Getting only CDS
  MEGA_CDS <- MEGA[MEGA$Gene_Type == "CDS",]
  # Getting only Synonymous
  MEGA_CDS_Syn <- MEGA_CDS[MEGA_CDS$Synonymous == "Synonymous",]
  # Subsetting the GTF
  V170_REF_GTF_CDS <- V170_REF_GTF[V170_REF_GTF$Gene_Type == "CDS",]
  # Function to get gene count !!! Do not use mut freq here sinc it has to be computed using numner of synonymous sites
  # See script R_Count_S_NS_Sites.R for more info
  mut_freq_gene(MEGA_CDS_Syn,V170_REF_GTF_CDS)
  # Load the counts for NSyn and Syn sites per gene
  Syn_sites_genes <- read.table("R_txt_files/Syn_NSyn_Sites_count.txt",header = T,sep="\t")
  # Using it to get the correct frequencies
  Mut_freqs_MEGA_CDS_Syn$MutationFreq <- Mut_freqs_MEGA_CDS_Syn$MutationCount / Syn_sites_genes$Syn
  Mut_freqs_MEGA_CDS_Syn$SIZE <- Syn_sites_genes$Syn
  colnames(Mut_freqs_MEGA_CDS_Syn) <- c("Gene_ID","MutationCount","START","Nb_Syn_Sites","Gene_Type","MutationFreq" )
}

##### For CDS NSyn #####
{   
  # Getting only CDS
  MEGA_CDS <- MEGA[MEGA$Gene_Type == "CDS",]
  # Getting only Synonymous
  MEGA_CDS_NSyn <- MEGA_CDS[MEGA_CDS$Synonymous == "NON_Synonymous",]
  # Subsetting the GTF
  V170_REF_GTF_CDS <- V170_REF_GTF[V170_REF_GTF$Gene_Type == "CDS",]
  # Function to get gene count !!! Do not use mut freq here sinc it has to be computed using numner of synonymous sites
  # See script R_Count_S_NS_Sites.R for more info
  mut_freq_gene(MEGA_CDS_NSyn,V170_REF_GTF_CDS)
  # Load the counts for NSyn and Syn sites per gene
  Syn_sites_genes <- read.table("R_txt_files/Syn_NSyn_Sites_count.txt",header = T,sep="\t")
  # Using it to get the correct frequencies
  Mut_freqs_MEGA_CDS_NSyn$MutationFreq <- Mut_freqs_MEGA_CDS_NSyn$MutationCount / Syn_sites_genes$NSyn
  Mut_freqs_MEGA_CDS_NSyn$SIZE <- Syn_sites_genes$NSyn
  colnames(Mut_freqs_MEGA_CDS_NSyn) <- c("Gene_ID","MutationCount","START","Nb_Syn_Sites","Gene_Type","MutationFreq" )
  }

##### Grouping into ETC complexes #####
{
  # Load the GTF
  V170_CDS <- V170_REF_GTF[V170_REF_GTF$Gene_Type == "CDS",]
  V170_CDS$SIZE <- V170_CDS$END - V170_CDS$START + 1
  
  # List of genes per ETC
  {
    ETC1 = c("ND6","ND4L","ND1","ND2","ND4","ND3","ND5")
    ETC3 = "CTB1"
    ETC4 = c("COX3","COX1","COX2")
    ETC5 = "ATP6"
    ALL = unique(V170_CDS$Gene_ID)
  }

  # Subset the mutation table
  {
    MEGA_ETC1 = MEGA[MEGA$Gene_ID %in% ETC1,]
    MEGA_ETC3 = MEGA[MEGA$Gene_ID %in% ETC3,]
    MEGA_ETC4 = MEGA[MEGA$Gene_ID %in% ETC4,]
    MEGA_ETC5 = MEGA[MEGA$Gene_ID %in% ETC5,]
    MEGA_ALL = MEGA[MEGA$Gene_ID %in% ALL,]
  }

  # Calculate the size of each complex
  {  
    ## Load GTF and name the columns
    V170_REF_GTF <- read.table("../01-DATA/9-VARIANT_CALLING/GTF/GTF_V170.fa.gtf")
    colnames(V170_REF_GTF) <- c("Chr","START","END","Gene_ID","Gene_Type")
    # Need to add a column with size
    V170_REF_GTF$SIZE <- V170_REF_GTF$END - V170_REF_GTF$START + 1
    V170_CDS <- V170_REF_GTF[V170_REF_GTF$Gene_Type == "CDS",]
    size_ETC1 = sum(V170_CDS$SIZE[V170_CDS$Gene_ID %in% ETC1])
    size_ETC3 = sum(V170_CDS$SIZE[V170_CDS$Gene_ID %in% ETC3])
    size_ETC4 = sum(V170_CDS$SIZE[V170_CDS$Gene_ID %in% ETC4])
    size_ETC5 = sum(V170_CDS$SIZE[V170_CDS$Gene_ID %in% ETC5])
    size_ALL = sum(V170_CDS$SIZE)
  }
  
  # Calculating mutation frequencies
  {
    mut_freq_ETC1 = nrow(MEGA_ETC1) / size_ETC1
    mut_freq_ETC3 = nrow(MEGA_ETC3) / size_ETC3
    mut_freq_ETC4 = nrow(MEGA_ETC4) / size_ETC4
    mut_freq_ETC5 = nrow(MEGA_ETC5) / size_ETC5
    mut_freq_ALL = nrow(MEGA_ALL) / size_ALL
    
    sum_mut_freq_ETC = sum(mut_freq_ETC1,mut_freq_ETC3,mut_freq_ETC4,mut_freq_ETC5)
  }
  
  # Proportion of Syn Non Syn
  {
    table(MEGA_ETC1$Synonymous) / sum(table(MEGA_ETC1$Synonymous))
    table(MEGA_ETC3$Synonymous) / sum(table(MEGA_ETC3$Synonymous))
    table(MEGA_ETC4$Synonymous) / sum(table(MEGA_ETC4$Synonymous))
    table(MEGA_ETC5$Synonymous) / sum(table(MEGA_ETC5$Synonymous))
    table(MEGA_ALL$Synonymous) / sum(table(MEGA_ALL$Synonymous))
  }
  
  # Number of Syn and NSyn sites
  {
    # Load the Syn site 
    Syn_sites_genes <- read.table("R_txt_files/Syn_NSyn_Sites_count.txt",header = T,sep="\t")
    
    Syn_site_ETC1 = Syn_sites_genes[Syn_sites_genes$Gene %in% ETC1,]
    Syn_site_ETC3 = Syn_sites_genes[Syn_sites_genes$Gene %in% ETC3,]
    Syn_site_ETC4 = Syn_sites_genes[Syn_sites_genes$Gene %in% ETC4,]
    Syn_site_ETC5 = Syn_sites_genes[Syn_sites_genes$Gene %in% ETC5,]
    
  }

}

##### Statistical tests #####
{
  ### Load before test ###
  {
    ## Load mutation table
    MEGA <- read.table("R_txt_files/Mutation_List.txt",header=T,sep="\t")
    
    # List of genes per ETC
    {
      ETC1 = c("ND6","ND4L","ND1","ND2","ND4","ND3","ND5")
      ETC3 = "CTB1"
      ETC4 = c("COX3","COX1","COX2")
      ETC5 = "ATP6"
    }
    
    # Subset the mutation table
    {
      MEGA_ETC1 = MEGA[MEGA$Gene_ID %in% ETC1,]
      MEGA_ETC3 = MEGA[MEGA$Gene_ID %in% ETC3,]
      MEGA_ETC4 = MEGA[MEGA$Gene_ID %in% ETC4,]
      MEGA_ETC5 = MEGA[MEGA$Gene_ID %in% ETC5,]
    }
    
    # Calculate the size of each complex
    {  
      ## Load GTF and name the columns
      V170_REF_GTF <- read.table("../01-DATA/9-VARIANT_CALLING/GTF/GTF_V170.fa.gtf")
      colnames(V170_REF_GTF) <- c("Chr","START","END","Gene_ID","Gene_Type")
      # Need to add a column with size
      V170_REF_GTF$SIZE <- V170_REF_GTF$END - V170_REF_GTF$START + 1
      V170_CDS <- V170_REF_GTF[V170_REF_GTF$Gene_Type == "CDS",]
      size_ETC1 = sum(V170_CDS$SIZE[V170_CDS$Gene_ID %in% ETC1])
      size_ETC3 = sum(V170_CDS$SIZE[V170_CDS$Gene_ID %in% ETC3])
      size_ETC4 = sum(V170_CDS$SIZE[V170_CDS$Gene_ID %in% ETC4])
      size_ETC5 = sum(V170_CDS$SIZE[V170_CDS$Gene_ID %in% ETC5])
      size_ALL = sum(V170_CDS$SIZE)
    }
  }
  
  ##### Mutation freq per gene #####
  {
    # Only getting CDS muts
    MEGA_CDS <- MEGA[MEGA$Gene_Type == "CDS",]
    # Only getting CDS from the GTF
    V170_REF_GTF_CDS <- V170_REF_GTF[V170_REF_GTF$Gene_Type == "CDS",]
    ## Here I am doing g tests with observed against expected if homogeneous
    
    # For all mutations
    g_test_per_gene_plus_post(MEGA_CDS,V170_REF_GTF_CDS)
    
    # For Syn_mutations
    MEGA_CDS_Syn <- MEGA_CDS[MEGA_CDS$Synonymous == "Synonymous",]
    g_test_per_gene_plus_post(MEGA_CDS_Syn,V170_REF_GTF_CDS)
    
    # For NSyn mutations
    MEGA_CDS_NSyn <- MEGA_CDS[MEGA_CDS$Synonymous == "NON_Synonymous",]
    g_test_per_gene_plus_post(MEGA_CDS_NSyn,V170_REF_GTF_CDS)
    
    # Doing a proportion test for NS:S
    # I use previously computed expected for both
    {
      NSS.test <- function(Gene) {
        # Creating the matrix for the fisher test
        mat <- matrix(ncol=2,nrow=2)
        # First line has observed counts for Nsyn and Syn mutations
        mat[1,] <- c(Gtable_MEGA_CDS_NSyn[Gene,"Count"],Gtable_MEGA_CDS_Syn[Gene,"Count"])
        # Second line has the expected under homogeneity
        mat[2,] <- c(Gtable_MEGA_CDS_NSyn[Gene,"Expected_count"],Gtable_MEGA_CDS_Syn[Gene,"Expected_count"])
        # Printing the Gene we are looking at
        print(Gene)
        # Executing the fisher test
        print(fisher.test(mat))
      }
      
      # Now that the function is made I use it on all the genes
      for (G in Gtable_MEGA_CDS$Gene_ID){ NSS.test(G) }
    }
  }
  
  ##### Mutation freq per ETC complex #####
  {
    #### For all mutations ####
    {
      # Get mut counts
      Count_ETC1 <- nrow(MEGA_ETC1)
      Count_ETC3 <- nrow(MEGA_ETC3)
      Count_ETC4 <- nrow(MEGA_ETC4)
      Count_ETC5 <- nrow(MEGA_ETC5) 
      Count_total <- sum(Count_ETC1,Count_ETC3,Count_ETC4,Count_ETC5)
      
      # Compute overall mut freq
      mut_freq_allETC <- Count_total/size_ALL
      
      # Compute expected based on size * overall mut freq
      expected_ETC1 <- mut_freq_allETC * size_ETC1
      expected_ETC3 <- mut_freq_allETC * size_ETC3
      expected_ETC4 <- mut_freq_allETC * size_ETC4
      expected_ETC5 <- mut_freq_allETC * size_ETC5
      
      # Compute expected prop
      expected_prop_ETC1 <- expected_ETC1/Count_total
      expected_prop_ETC3 <- expected_ETC3/Count_total
      expected_prop_ETC4 <- expected_ETC4/Count_total
      expected_prop_ETC5 <- expected_ETC5/Count_total
      
      # Compute g-test of observed counts against expected proportions with Bonferoni correction
      G <- g.test(c(Count_ETC1,Count_ETC3,Count_ETC4,Count_ETC5),p=c(expected_prop_ETC1,expected_prop_ETC3,expected_prop_ETC4,expected_prop_ETC5))
      print(G)
      
      # Bonferroni correction = significance level / nb of comparison
      Bonferroni = 0.05/4
      
      # ETC specific test, compare each complex with the rest of the genome
      expected_prop_noETC1 <- 1-expected_prop_ETC1
      expected_prop_noETC3 <- 1-expected_prop_ETC3
      expected_prop_noETC4 <- 1-expected_prop_ETC4
      expected_prop_noETC5 <- 1-expected_prop_ETC5
      Count_noETC1 <- Count_total - Count_ETC1
      Count_noETC3 <- Count_total - Count_ETC3
      Count_noETC4 <- Count_total - Count_ETC4
      Count_noETC5 <- Count_total - Count_ETC5
      
      g.test(c(Count_ETC1,Count_noETC1),p=c(expected_prop_ETC1,expected_prop_noETC1))
      g.test(c(Count_ETC3,Count_noETC3),p=c(expected_prop_ETC3,expected_prop_noETC3))
      g.test(c(Count_ETC4,Count_noETC4),p=c(expected_prop_ETC4,expected_prop_noETC4))
      g.test(c(Count_ETC5,Count_noETC5),p=c(expected_prop_ETC5,expected_prop_noETC5))
    }
    
    #### For Syn_mutations (!!! Careful I use very similar script as before so run separaterly) ####
    {
      # Subset table
      MEGA_Syn_ETC1 <- MEGA_ETC1[MEGA_ETC1$Synonymous == "Synonymous",]
      MEGA_Syn_ETC3 <- MEGA_ETC3[MEGA_ETC3$Synonymous == "Synonymous",]
      MEGA_Syn_ETC4 <- MEGA_ETC4[MEGA_ETC4$Synonymous == "Synonymous",]
      MEGA_Syn_ETC5 <- MEGA_ETC5[MEGA_ETC5$Synonymous == "Synonymous",]
      
      # Get mut counts
      Count_ETC1 <- nrow(MEGA_Syn_ETC1)
      Count_ETC3 <- nrow(MEGA_Syn_ETC3)
      Count_ETC4 <- nrow(MEGA_Syn_ETC4)
      Count_ETC5 <- nrow(MEGA_Syn_ETC5) 
      Count_total <- sum(Count_ETC1,Count_ETC3,Count_ETC4,Count_ETC5)
      
      # Compute overall mut freq
      mut_freq_allETC <- Count_total/size_ALL
      
      # Compute expected based on size * overall mut freq
      expected_ETC1 <- mut_freq_allETC * size_ETC1
      expected_ETC3 <- mut_freq_allETC * size_ETC3
      expected_ETC4 <- mut_freq_allETC * size_ETC4
      expected_ETC5 <- mut_freq_allETC * size_ETC5
      
      # Compute expected prop
      expected_prop_ETC1 <- expected_ETC1/Count_total
      expected_prop_ETC3 <- expected_ETC3/Count_total
      expected_prop_ETC4 <- expected_ETC4/Count_total
      expected_prop_ETC5 <- expected_ETC5/Count_total
      
      # Compute g-test of observed counts against expected proportions with Bonferoni correction
      g.test(c(Count_ETC1,Count_ETC3,Count_ETC4,Count_ETC5),p=c(expected_prop_ETC1,expected_prop_ETC3,expected_prop_ETC4,expected_prop_ETC5))
      
      # Bonferroni correction = significance level / nb of comparison
      Bonferroni = 0.05/4
      
      # ETC specific test, compare each complex with the rest of the genome
      expected_prop_noETC1 <- 1-expected_prop_ETC1
      expected_prop_noETC3 <- 1-expected_prop_ETC3
      expected_prop_noETC4 <- 1-expected_prop_ETC4
      expected_prop_noETC5 <- 1-expected_prop_ETC5
      Count_noETC1 <- Count_total - Count_ETC1
      Count_noETC3 <- Count_total - Count_ETC3
      Count_noETC4 <- Count_total - Count_ETC4
      Count_noETC5 <- Count_total - Count_ETC5
      
      g.test(c(Count_ETC1,Count_noETC1),p=c(expected_prop_ETC1,expected_prop_noETC1))
      g.test(c(Count_ETC3,Count_noETC3),p=c(expected_prop_ETC3,expected_prop_noETC3))
      g.test(c(Count_ETC4,Count_noETC4),p=c(expected_prop_ETC4,expected_prop_noETC4))
      g.test(c(Count_ETC5,Count_noETC5),p=c(expected_prop_ETC5,expected_prop_noETC5))
    }
    
    #### For NSyn mutations (!!! Careful I use very similar script as before so run separaterly) ####
    {
      
      # Subset table
      MEGA_NSyn_ETC1 <- MEGA_ETC1[MEGA_ETC1$Synonymous == "NON_Synonymous",]
      MEGA_NSyn_ETC3 <- MEGA_ETC3[MEGA_ETC3$Synonymous == "NON_Synonymous",]
      MEGA_NSyn_ETC4 <- MEGA_ETC4[MEGA_ETC4$Synonymous == "NON_Synonymous",]
      MEGA_NSyn_ETC5 <- MEGA_ETC5[MEGA_ETC5$Synonymous == "NON_Synonymous",]
      
      # Get mut counts
      Count_ETC1 <- nrow(MEGA_NSyn_ETC1)
      Count_ETC3 <- nrow(MEGA_NSyn_ETC3)
      Count_ETC4 <- nrow(MEGA_NSyn_ETC4)
      Count_ETC5 <- nrow(MEGA_NSyn_ETC5) 
      Count_total <- sum(Count_ETC1,Count_ETC3,Count_ETC4,Count_ETC5)
      
      # Compute overall mut freq
      mut_freq_allETC <- Count_total/size_ALL
      
      # Compute expected based on size * overall mut freq
      expected_ETC1 <- mut_freq_allETC * size_ETC1
      expected_ETC3 <- mut_freq_allETC * size_ETC3
      expected_ETC4 <- mut_freq_allETC * size_ETC4
      expected_ETC5 <- mut_freq_allETC * size_ETC5
      
      # Compute expected prop
      expected_prop_ETC1 <- expected_ETC1/Count_total
      expected_prop_ETC3 <- expected_ETC3/Count_total
      expected_prop_ETC4 <- expected_ETC4/Count_total
      expected_prop_ETC5 <- expected_ETC5/Count_total
      
      # Compute g-test of observed counts against expected proportions with Bonferoni correction
      g.test(c(Count_ETC1,Count_ETC3,Count_ETC4,Count_ETC5),p=c(expected_prop_ETC1,expected_prop_ETC3,expected_prop_ETC4,expected_prop_ETC5))
      
      # Bonferroni correction = significance level / nb of comparison
      Bonferroni = 0.05/4
      
      # ETC specific test, compare each complex with the rest of the genome
      expected_prop_noETC1 <- 1-expected_prop_ETC1
      expected_prop_noETC3 <- 1-expected_prop_ETC3
      expected_prop_noETC4 <- 1-expected_prop_ETC4
      expected_prop_noETC5 <- 1-expected_prop_ETC5
      Count_noETC1 <- Count_total - Count_ETC1
      Count_noETC3 <- Count_total - Count_ETC3
      Count_noETC4 <- Count_total - Count_ETC4
      Count_noETC5 <- Count_total - Count_ETC5
      
      g.test(c(Count_ETC1,Count_noETC1),p=c(expected_prop_ETC1,expected_prop_noETC1))
      g.test(c(Count_ETC3,Count_noETC3),p=c(expected_prop_ETC3,expected_prop_noETC3))
      g.test(c(Count_ETC4,Count_noETC4),p=c(expected_prop_ETC4,expected_prop_noETC4))
      g.test(c(Count_ETC5,Count_noETC5),p=c(expected_prop_ETC5,expected_prop_noETC5))
    }
    
    #### For Syn_mutations WITHOUT ETC1 (!!! Careful I use very similar script as before so run separaterly) ####
    {
      # Subset table
      MEGA_Syn_ETC1 <- MEGA_ETC1[MEGA_ETC1$Synonymous == "Synonymous",]
      MEGA_Syn_ETC3 <- MEGA_ETC3[MEGA_ETC3$Synonymous == "Synonymous",]
      MEGA_Syn_ETC4 <- MEGA_ETC4[MEGA_ETC4$Synonymous == "Synonymous",]
      MEGA_Syn_ETC5 <- MEGA_ETC5[MEGA_ETC5$Synonymous == "Synonymous",]
      
      # Get mut counts
      Count_ETC1 <- nrow(MEGA_Syn_ETC1)
      Count_ETC3 <- nrow(MEGA_Syn_ETC3)
      Count_ETC4 <- nrow(MEGA_Syn_ETC4)
      Count_ETC5 <- nrow(MEGA_Syn_ETC5) 
      Count_total <- sum(Count_ETC3,Count_ETC4,Count_ETC5)
      
      # Compute overall mut freq
      mut_freq_allETC <- Count_total/(size_ALL - size_ETC1)
      
      # Compute expected based on size * overall mut freq
      expected_ETC1 <- mut_freq_allETC * size_ETC1
      expected_ETC3 <- mut_freq_allETC * size_ETC3
      expected_ETC4 <- mut_freq_allETC * size_ETC4
      expected_ETC5 <- mut_freq_allETC * size_ETC5
      
      # Compute expected prop
      expected_prop_ETC1 <- expected_ETC1/Count_total
      expected_prop_ETC3 <- expected_ETC3/Count_total
      expected_prop_ETC4 <- expected_ETC4/Count_total
      expected_prop_ETC5 <- expected_ETC5/Count_total
      
      # Compute g-test of observed counts against expected proportions with Bonferoni correction
      g.test(c(Count_ETC3,Count_ETC4,Count_ETC5),p=c(expected_prop_ETC3,expected_prop_ETC4,expected_prop_ETC5))
      
      # Bonferroni correction = significance level / nb of comparison
      Bonferroni = 0.05/4
      
      # Fix expected prop just for ETC1
      expected_prop_ETC1 <- expected_prop_ETC1/sum(expected_prop_ETC1,expected_prop_ETC3,expected_prop_ETC4,expected_prop_ETC5)
      
      # ETC specific test, compare each complex with the rest of the genome
      expected_prop_noETC1 <- 1-expected_prop_ETC1
      expected_prop_noETC3 <- 1-expected_prop_ETC3
      expected_prop_noETC4 <- 1-expected_prop_ETC4
      expected_prop_noETC5 <- 1-expected_prop_ETC5
      
      Count_noETC3 <- Count_total - Count_ETC3
      Count_noETC4 <- Count_total - Count_ETC4
      Count_noETC5 <- Count_total - Count_ETC5
      Count_noETC1 <- sum(Count_ETC3,Count_ETC4,Count_ETC5)
      

      g.test(c(Count_ETC1,Count_noETC1),p=c(expected_prop_ETC1,expected_prop_noETC1))
      g.test(c(Count_ETC3,Count_noETC3),p=c(expected_prop_ETC3,expected_prop_noETC3))
      g.test(c(Count_ETC4,Count_noETC4),p=c(expected_prop_ETC4,expected_prop_noETC4))
      g.test(c(Count_ETC5,Count_noETC5),p=c(expected_prop_ETC5,expected_prop_noETC5))
    }
    
    #### For NSyn mutations WITHOUT ETC1 (!!! Careful I use very similar script as before so run separaterly) ####
    {
      # Subset table
      MEGA_NSyn_ETC1 <- MEGA_ETC1[MEGA_ETC1$Synonymous == "NON_Synonymous",]
      MEGA_NSyn_ETC3 <- MEGA_ETC3[MEGA_ETC3$Synonymous == "NON_Synonymous",]
      MEGA_NSyn_ETC4 <- MEGA_ETC4[MEGA_ETC4$Synonymous == "NON_Synonymous",]
      MEGA_NSyn_ETC5 <- MEGA_ETC5[MEGA_ETC5$Synonymous == "NON_Synonymous",]
      
      # Get mut counts
      Count_ETC1 <- nrow(MEGA_NSyn_ETC1)
      Count_ETC3 <- nrow(MEGA_NSyn_ETC3)
      Count_ETC4 <- nrow(MEGA_NSyn_ETC4)
      Count_ETC5 <- nrow(MEGA_NSyn_ETC5) 
      Count_total <- sum(Count_ETC3,Count_ETC4,Count_ETC5)
      
      # Compute overall mut freq
      mut_freq_allETC <- Count_total/(size_ALL - size_ETC1)
      
      # Compute expected based on size * overall mut freq
      expected_ETC1 <- mut_freq_allETC * size_ETC1
      expected_ETC3 <- mut_freq_allETC * size_ETC3
      expected_ETC4 <- mut_freq_allETC * size_ETC4
      expected_ETC5 <- mut_freq_allETC * size_ETC5
      
      # Compute expected prop
      expected_prop_ETC1 <- expected_ETC1/Count_total
      expected_prop_ETC3 <- expected_ETC3/Count_total
      expected_prop_ETC4 <- expected_ETC4/Count_total
      expected_prop_ETC5 <- expected_ETC5/Count_total
      
      # Compute g-test of observed counts against expected proportions with Bonferoni correction
      g.test(c(Count_ETC3,Count_ETC4,Count_ETC5),p=c(expected_prop_ETC3,expected_prop_ETC4,expected_prop_ETC5))
      
      # Bonferroni correction = significance level / nb of comparison
      Bonferroni = 0.05/4
      
      # Fix expected prop just for ETC1
      expected_prop_ETC1 <- expected_prop_ETC1/sum(expected_prop_ETC1,expected_prop_ETC3,expected_prop_ETC4,expected_prop_ETC5)
      
      # ETC specific test, compare each complex with the rest of the genome
      expected_prop_noETC1 <- 1-expected_prop_ETC1
      expected_prop_noETC3 <- 1-expected_prop_ETC3
      expected_prop_noETC4 <- 1-expected_prop_ETC4
      expected_prop_noETC5 <- 1-expected_prop_ETC5
      
      Count_noETC3 <- Count_total - Count_ETC3
      Count_noETC4 <- Count_total - Count_ETC4
      Count_noETC5 <- Count_total - Count_ETC5
      Count_noETC1 <- sum(Count_ETC3,Count_ETC4,Count_ETC5)
      
      
      g.test(c(Count_ETC1,Count_noETC1),p=c(expected_prop_ETC1,expected_prop_noETC1))
      g.test(c(Count_ETC3,Count_noETC3),p=c(expected_prop_ETC3,expected_prop_noETC3))
      g.test(c(Count_ETC4,Count_noETC4),p=c(expected_prop_ETC4,expected_prop_noETC4))
      g.test(c(Count_ETC5,Count_noETC5),p=c(expected_prop_ETC5,expected_prop_noETC5))
    }
    
    #### Test for NS:S fisher.test comparing observed and expected counts ####
    {
      ## Run the g.test.per.gene for NSyn and Syn mutations above
      # Getting the expected counts for NSyn
      expected_NSyn_ETC1 <- sum(Gtable_MEGA_CDS_NSyn$Expected_count[Gtable_MEGA_CDS_NSyn$Gene_ID %in% ETC1])
      expected_NSyn_ETC3 <- sum(Gtable_MEGA_CDS_NSyn$Expected_count[Gtable_MEGA_CDS_NSyn$Gene_ID %in% ETC3])
      expected_NSyn_ETC4 <- sum(Gtable_MEGA_CDS_NSyn$Expected_count[Gtable_MEGA_CDS_NSyn$Gene_ID %in% ETC4])
      expected_NSyn_ETC5 <- sum(Gtable_MEGA_CDS_NSyn$Expected_count[Gtable_MEGA_CDS_NSyn$Gene_ID %in% ETC5])
      # Getting the expected counts for Syn
      expected_Syn_ETC1 <- sum(Gtable_MEGA_CDS_Syn$Expected_count[Gtable_MEGA_CDS_Syn$Gene_ID %in% ETC1])
      expected_Syn_ETC3 <- sum(Gtable_MEGA_CDS_Syn$Expected_count[Gtable_MEGA_CDS_Syn$Gene_ID %in% ETC3])
      expected_Syn_ETC4 <- sum(Gtable_MEGA_CDS_Syn$Expected_count[Gtable_MEGA_CDS_Syn$Gene_ID %in% ETC4])
      expected_Syn_ETC5 <- sum(Gtable_MEGA_CDS_Syn$Expected_count[Gtable_MEGA_CDS_Syn$Gene_ID %in% ETC5])
      # Getting the observed count for NSyn
      observed_NSyn_ETC1 <- sum(Gtable_MEGA_CDS_NSyn$Count[Gtable_MEGA_CDS_NSyn$Gene_ID %in% ETC1])
      observed_NSyn_ETC3 <- sum(Gtable_MEGA_CDS_NSyn$Count[Gtable_MEGA_CDS_NSyn$Gene_ID %in% ETC3])
      observed_NSyn_ETC4 <- sum(Gtable_MEGA_CDS_NSyn$Count[Gtable_MEGA_CDS_NSyn$Gene_ID %in% ETC4])
      observed_NSyn_ETC5 <- sum(Gtable_MEGA_CDS_NSyn$Count[Gtable_MEGA_CDS_NSyn$Gene_ID %in% ETC5])
      # Getting the observed count for Syn
      observed_Syn_ETC1 <- sum(Gtable_MEGA_CDS_Syn$Count[Gtable_MEGA_CDS_Syn$Gene_ID %in% ETC1])
      observed_Syn_ETC3 <- sum(Gtable_MEGA_CDS_Syn$Count[Gtable_MEGA_CDS_Syn$Gene_ID %in% ETC3])
      observed_Syn_ETC4 <- sum(Gtable_MEGA_CDS_Syn$Count[Gtable_MEGA_CDS_Syn$Gene_ID %in% ETC4])
      observed_Syn_ETC5 <- sum(Gtable_MEGA_CDS_Syn$Count[Gtable_MEGA_CDS_Syn$Gene_ID %in% ETC5])
      # Building the matrix
      mat_ETC1 <- matrix(c(observed_NSyn_ETC1,expected_NSyn_ETC1,observed_Syn_ETC1,expected_Syn_ETC1),nrow=2,ncol=2)
      mat_ETC3 <- matrix(c(observed_NSyn_ETC3,expected_NSyn_ETC3,observed_Syn_ETC3,expected_Syn_ETC3),nrow=2,ncol=2)
      mat_ETC4 <- matrix(c(observed_NSyn_ETC4,expected_NSyn_ETC4,observed_Syn_ETC4,expected_Syn_ETC4),nrow=2,ncol=2)
      mat_ETC5 <- matrix(c(observed_NSyn_ETC5,expected_NSyn_ETC5,observed_Syn_ETC5,expected_Syn_ETC5),nrow=2,ncol=2)
      # Fisher test
      fisher.test(mat_ETC1)
      fisher.test(mat_ETC3)
      fisher.test(mat_ETC4)
      fisher.test(mat_ETC5)
      }
    
  }
}









