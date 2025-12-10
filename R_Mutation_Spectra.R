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
}

###### Load the GTF I use ######
{
  V170_REF_GTF = read.table("../01-DATA/9-VARIANT_CALLING/GTF/GTF_V170.fa.gtf")
  V170_REF_GTF <- subset(V170_REF_GTF, select = -V1)
  V170_REF_GTF[,5] <- V170_REF_GTF[,2] - V170_REF_GTF[,1] +1
  colnames(V170_REF_GTF) = c("START","END","Gene_ID","Gene_Type","SIZE")
}

###### Plotting Mutation spectra ######
{
  # Load the reference fasta (needed to count number of bp)
  library(seqinr)
  V170 <- read.fasta("../01-DATA/9-VARIANT_CALLING/REF/V170.fa")
  
  ###### Plot for all mutations (look at the table Spectrum_MEGA) #####
  spectrum_plot(MEGA,REF = V170[[1]])
  
  ###### Plot for four-fold degenerate sites  ###### 
  {
    # I have to use base composition at 4x sites so Corrected values here cannot be used
    MEGA_4x <- read.table("R_txt_files/Mutation_Fourfold.txt",sep="\t",header=T)
    spectrum_plot(MEGA_4x, REF=V170[[1]],Corrected_Bases = F)
  }
  
  ###### Plot for two-fold degenerate sites  ###### 
  {
    # I have to use base composition at 4x sites so Corrected values here cannot be used
    MEGA_2x <- read.table("R_txt_files/Mutation_Twofold.txt",sep="\t",header=T)
    spectrum_plot(MEGA_2x, REF=V170[[1]],Corrected_Bases = F)
  }
  
  ###### Plot for non-degenerate sites  ###### 
  {
    # I have to use base composition at 4x sites so Corrected values here cannot be used
    MEGA_0x <- read.table("R_txt_files/Mutation_NonDegenerate.txt",sep="\t",header=T)
    spectrum_plot(MEGA_0x, REF=V170[[1]],Corrected_Bases = F)
  }

  ##### Plot for Leuthner data #####
  {
    # Because the formating of their data sucks I have to do it that way and cannot use my function
    mut_name <- c("A -> T","A -> C","A -> G","T -> A","T -> C","T -> G","C -> A","C -> T","C -> G","G -> A","G -> T","G -> C")
    mut_counts <- c(28,7,19,48,25,15,102,37,98,56,198,127)
    Leuthner_temp = data.frame(Mutation = mut_name, Count = mut_counts)
    row.names(Leuthner_temp) <- mut_name
    
    # Creating the proper table
    list.spectrum = c("A:T -> C:G","A:T -> T:A","C:G -> A:T","C:G -> G:C","A:T -> G:C","C:G -> T:A")
    list.transi = c("A:T -> G:C","C:G -> T:A")
    Leuthner = data.frame(Mutation=list.spectrum,Count=0)
    row.names(Leuthner) <- list.spectrum
    
    # Grouping into the mut spectrum mutations
    Leuthner["A:T -> C:G","Count"] <- sum(Leuthner_temp["A -> C","Count"],Leuthner_temp["T -> G","Count"])
    Leuthner["A:T -> T:A","Count"] <- sum(Leuthner_temp["A -> T","Count"],Leuthner_temp["T -> A","Count"])
    Leuthner["C:G -> A:T","Count"] <- sum(Leuthner_temp["G -> T","Count"],Leuthner_temp["C -> A","Count"])
    Leuthner["C:G -> G:C","Count"] <- sum(Leuthner_temp["G -> C","Count"],Leuthner_temp["C -> G","Count"])
    Leuthner["A:T -> G:C","Count"] <- sum(Leuthner_temp["A -> G","Count"],Leuthner_temp["T -> C","Count"])
    Leuthner["C:G -> T:A","Count"] <- sum(Leuthner_temp["G -> A","Count"],Leuthner_temp["C -> T","Count"])
    
    # Removing Leuthner_temp
    rm(Leuthner_temp)
    
    # Adding transition info
    Leuthner$Transition <- ""
    Leuthner$Transition[Leuthner$Mutation %in% list.transi] <- "Transition"
    Leuthner$Transition[Leuthner$Mutation %in% list.transi == F] <- "Transversion"
    
    ## Calculating proportions
    Leuthner$Proportion <- Leuthner$Count / sum(Leuthner$Count)
    
    ## Counting base composition to account for it
    library(seqinr)
    V170 <- read.fasta("../01-DATA/9-VARIANT_CALLING/REF/V170.fa")
    Bases_REF = seqinr::count(V170[[1]],wordsize = 1)
    
    ## Correct by number of bases
    nbGC = Bases_REF[["g"]] + Bases_REF[["c"]]
    nbAT = Bases_REF[["a"]] + Bases_REF[["t"]]
    
    Leuthner[c("A:T -> C:G","A:T -> T:A","A:T -> G:C"),"Corrected"] <- Leuthner[c("A:T -> C:G","A:T -> T:A","A:T -> G:C"),"Count"] / nbAT
    Leuthner[c("C:G -> A:T","C:G -> G:C","C:G -> T:A"),"Corrected"] <- Leuthner[c("C:G -> A:T","C:G -> G:C","C:G -> T:A"),"Count"] / nbGC
    
    ## Percentage
    Leuthner$Normalized <- Leuthner$Corrected * (1/sum(Leuthner$Corrected))
    
    ## ggplot
    library(ggplot2)
    ggplot(Leuthner, aes(x=Mutation, y=Normalized, fill=Transition)) +
      geom_col(colour="black") + theme_minimal() + labs(title = paste("Mutation spectrum for",deparse(substitute(Table))))+ 
      theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5), legend.title = element_blank())+ 
      scale_fill_manual(values = c("cadetblue","grey90"))+
      ylim(0,1)+
      ylab("Percentage")
  }
  
  ##### Plot for Waneka data #####
  {
    library(openxlsx)
    Waneka = read.xlsx("Previous_results/Waneka/File_S1.xlsx")
    library(stringr)
    Waneka$Mutation <- str_replace(Waneka$`single.stranded.substitution.type.on.F-strand`, ">", " -> ")
    spectrum_plot(Waneka)
  }
  
  ##### Plot for Konrad's data #####
  {
    Konrad = read.table("Previous_results/Konrad/Table1_Konrad2017_MALines_raw.txt",header = T,sep="\t")
    # Only taking SNMs
    Konrad_SNM <- Konrad[c(1,2,3,5,7,12,22,23,24),]
  }
  
}

##### Statistical tests #####
{
  ##### Pairwise proportion test on observed mutations #####
  # The following vectors are the counts for each mutation type in the following order
  # A:T -> C:G | A:T -> T:A | C:G -> A:T | C:G -> G:C | A:T -> G:C | C:G -> T:A
  This <- c(54,224,144,9,1057,875)
  ThisFourfold <- c(42,133,53,5,408,312)
  Leuthner <- c(22,76,300,225,44,93)
  Waneka <- c(2,20,79,8,50,94)
  
  test_mut_counts <- function(x,y) {
    matrix <- rbind(x,y)
    library(AMR)
    print(g.test(matrix))
    print(fisher.test(matrix,simulate.p.value=T))
  }

  test_mut_counts(This,ThisFourfold)
  test_mut_counts(This,Waneka)
  test_mut_counts(This,Leuthner)
  test_mut_counts(ThisFourfold,Waneka)
  test_mut_counts(ThisFourfold,Leuthner)
  test_mut_counts(Leuthner,Waneka)
}