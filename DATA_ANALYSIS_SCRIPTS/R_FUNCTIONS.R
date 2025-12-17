### Function to make a spectrum plot out of a table of mutations ###
spectrum_plot <- function(Table,REF = V170[[1]], Corrected_Bases = T) {
  
  ## Making the spectrum table
  spectrum = data.frame((matrix(nrow=6, ncol=6)))
  colnames(spectrum) = c("Mutation","Count","Proportion","Corrected","Normalized","Transition")
  list.spectrum = c("A:T -> C:G","A:T -> T:A","C:G -> A:T","C:G -> G:C","A:T -> G:C","C:G -> T:A")
  list.transi = c("A:T -> G:C","C:G -> T:A")
  spectrum$Mutation = list.spectrum
  spectrum$Transition <- "Transversion"
  spectrum[which(spectrum$Mutation %in% list.transi),"Transition"] <- "Transition"
  row.names(spectrum) = list.spectrum
  
  ## ==> I count each type of mutation
  if ("MUTATION" %in% colnames(Table)) {
    spectrum["A:T -> C:G","Count"] <- sum(Table$MUTATION == "A -> C",Table$MUTATION == "T -> G")
    spectrum["A:T -> T:A","Count"] <- sum(Table$MUTATION == "A -> T",Table$MUTATION == "T -> A")
    spectrum["C:G -> A:T","Count"] <- sum(Table$MUTATION == "C -> A",Table$MUTATION == "G -> T")
    spectrum["C:G -> G:C","Count"] <- sum(Table$MUTATION == "C -> G",Table$MUTATION == "G -> C")
    spectrum["A:T -> G:C","Count"] <- sum(Table$MUTATION == "A -> G",Table$MUTATION == "T -> C")
    spectrum["C:G -> T:A","Count"] <- sum(Table$MUTATION == "C -> T",Table$MUTATION == "G -> A")

  } else {
    if ("Mutation" %in% colnames(Table)) {
      spectrum["A:T -> C:G","Count"] <- sum(Table$Mutation == "A -> C",Table$Mutation == "T -> G")
      spectrum["A:T -> T:A","Count"] <- sum(Table$Mutation == "A -> T",Table$Mutation == "T -> A")
      spectrum["C:G -> A:T","Count"] <- sum(Table$Mutation == "C -> A",Table$Mutation == "G -> T")
      spectrum["C:G -> G:C","Count"] <- sum(Table$Mutation == "C -> G",Table$Mutation == "G -> C")
      spectrum["A:T -> G:C","Count"] <- sum(Table$Mutation == "A -> G",Table$Mutation == "T -> C")
      spectrum["C:G -> T:A","Count"] <- sum(Table$Mutation == "C -> T",Table$Mutation == "G -> A")
    } else { print("Missing a collumn with SNM, called Mutation or MUTATION with values as: G -> A") }
  }
  
  spectrum$Proportion <- spectrum$Count / sum(spectrum$Count)
  
  ## Counting bases in a REF seq
  library(seqinr)
  Bases_REF = seqinr::count(REF,wordsize = 1)
  
  if (Corrected_Bases) {
    ## Correct by number of bases
    nbGC = Bases_REF[["g"]] + Bases_REF[["c"]]
    nbAT = Bases_REF[["a"]] + Bases_REF[["t"]]
    
    spectrum[c("A:T -> C:G","A:T -> T:A","A:T -> G:C"),"Corrected"] <- spectrum[c("A:T -> C:G","A:T -> T:A","A:T -> G:C"),"Count"] / nbAT
    spectrum[c("C:G -> A:T","C:G -> G:C","C:G -> T:A"),"Corrected"] <- spectrum[c("C:G -> A:T","C:G -> G:C","C:G -> T:A"),"Count"] / nbGC
    
    ## Percentage
    spectrum$Normalized <- spectrum$Corrected * (1/sum(spectrum$Corrected))
  }

  ## Print the table
  print(spectrum)
  print("Spectrum table can be found as global variable named Spectrun_[Table]")
  assign(paste("Spectrum_",deparse(substitute(Table)),sep=""),spectrum, envir = globalenv())
  
  if (Corrected_Bases) {
    ## ggplot
    library(ggplot2)
    ggplot(spectrum, aes(x=Mutation, y=Normalized, fill=Transition)) +
      geom_col(colour="black") + theme_minimal() + labs(title = paste("Mutation spectrum for",deparse(substitute(Table))))+ 
      theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5), legend.title = element_blank())+ 
      scale_fill_manual(values = c("cadetblue","grey90"))+
      ylim(0,1)+
      ylab("Percentage")
    }

}
####################################################################################################################

####################################################################################################################
### Function for mut rate per gene
# colnames(GTF) = c("START","END","Gene_ID","Gene_Type","SIZE")
mut_freq_gene <- function(Table,GTF = V170_REF_GTF) {
  
  library(dplyr)
  
  # Mutation counts 
  mutation_counts <- Table %>%
    group_by(Gene_ID) %>%
    summarize(MutationCount = n())
  
  # Create a complete list of Gene_IDs from V152_REF_GTF
  complete_gene_list <- unique(GTF$Gene_ID)
  
  # Merge the complete gene list with mutation counts and START
  merged_data <- data.frame(Gene_ID = complete_gene_list) %>%
    left_join(mutation_counts, by = "Gene_ID") %>%
    left_join(select(GTF, Gene_ID, START, SIZE, Gene_Type), by = "Gene_ID") %>%
    left_join(select(GTF, Gene_ID), by = "Gene_ID")
  
  # Replace missing MutationCount values with 0
  merged_data$MutationCount[is.na(merged_data$MutationCount)] <- 0
  
  # Calculate the mutation rate by dividing the mutation count by gene size
  merged_data$MutationFreq <- merged_data$MutationCount / merged_data$SIZE
  
  # mutation rate across the entire genome
  mut_rate_tot = nrow(Table) / sum(GTF$SIZE) 
  
  # Create a customized and visually appealing bar plot for mutation rate
  library(ggplot2)
  
  rate_line <- data.frame(yintercept = mut_rate_tot, Line = "Total mut freq")
  assign(paste("Mut_freqs_",deparse(substitute(Table)),sep=""),merged_data, envir = globalenv())
  print("Mutation freq overall")
  print(rate_line)
  print("#########################################")
  print(merged_data)
  
  ggplot(merged_data, aes(x = reorder(Gene_ID, START), y = MutationFreq, fill = Gene_Type)) +
    geom_bar(stat = "identity", color = "black") +
    labs(title = paste("Mutation frequency per gene for",deparse(substitute(Table))), x = "Gene ID", y = "Frequency (Count / Nb bp)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    geom_hline(aes(yintercept=yintercept, linetype=Line),rate_line,col="red")+
    scale_fill_brewer()
}
####################################################################################################################

####################################################################################################################
### Function g test from a mutation list and a GTF (specific to a set of genes)
# The list of genes in the GTF set the test comparisons
g_test_per_gene_plus_post <- function(Table,GTF) {
  
  library(AMR) 
  library(dplyr)
  
  ## Table of observed and expected for the GeneType provided
  genetype_size <- sum(GTF$SIZE)
  mut_freq <- nrow(Table) / genetype_size
  count = data.frame(table(Table$Gene_ID))
  gt_table <- data.frame(Gene_ID = count$Var1, Count = count$Freq) %>%
    right_join(GTF, by="Gene_ID")
  gt_table$Expected_count <- mut_freq*gt_table$SIZE
  gt_table$Expected_percentage <- gt_table$Expected_count / sum(gt_table$Expected_count)
  gt_table[is.na(gt_table)] <- 0
  
  G = g.test(gt_table$Count, p = gt_table$Expected_percentage)
  print("p-value for overall fit of our data distribution compared to expected if gene mut rate was equal")
  print(G$p.value)
  print(G$statistic)
  if (G$p.value > 0.05) { print("No differences from expected if gene mut rates were equal") } else {
    print("Significant difference between our distribution and what would be expected if mut rates were constent across genes")
    print("Next step is to look at each gene to see if it has a different mut rate compared to the rest of the genome")
    ## http://www.biostathandbook.com/exactgof.html #posthoc
    
    ## Bonferroni correction = significance level / nb of comparison
    Bonferroni = 0.05/length(unique(gt_table$Gene_ID))
    
    ## Add 3 columns
    gt_table$PVAL <- 0
    gt_table$Significative <- 0
    gt_table$Statistic <- 0
    gt_table$df <- 0
    gt_table$Bonferroni <- Bonferroni
    
    # rename rows as gene for easier access
    row.names(gt_table) <- gt_table$Gene_ID
    
    ## Loop for G test for each gene vs rest
    for (gene in gt_table$Gene_ID) {
      
      nb_mut_gene = gt_table[gene,"Count"]
      nb_mut_all_but = nrow(Table[Table$Gene_ID != gene,])
      obs_percentage_gene = gt_table$Expected_percentage[gt_table$Gene_ID == gene]
      exp_percentage_all_but = sum(gt_table$Expected_percentage[gt_table$Gene_ID != gene])
      
      nb <- c(nb_mut_gene,nb_mut_all_but)
      proba <- c(obs_percentage_gene,exp_percentage_all_but)
      
      G2 = g.test(nb,p=proba)
      gt_table[gene,"df"] <- G2[["argument"]][["df"]]
      gt_table[gene,"PVAL"] <- G2$p.value
      gt_table[gene,"Statistic"] <- G2$statistic
      if (gt_table[gene,"PVAL"] < Bonferroni) { 
        gt_table[gene,"Significative"] <- "Significative" } else { 
          gt_table[gene,"Significative"] <- "NOT" }
      
    }
    print(gt_table)
    print("Finale table is available as displayed in")
    print(paste("Gtable_",deparse(substitute(Table)),sep=""))
    assign(paste("Gtable_",deparse(substitute(Table)),sep=""), gt_table, envir = globalenv())  
  }
  
}
####################################################################################################################

####################################################################################################################
### Function for taging protein coding gene sequences positions by the degeneracy of their positions
# Input 1 sequence as a fasta, import the sequence with read.fasta (seqinr)
# The sequence as to be in a list or as a vector
# Can also be a list of vectors with multiple genes
# example with the ref sequence for c elegans CDS
# option "codon" to only do it for only one pos of a codon, codon can take F for False, 1, 2 or 3 for various codon position, the sequence given has to be a 3 character sequence (a codon)
#### !!!!!!!!!!!!!!!!!!   VERY IMPORTANT TO HAVE THE SEQ IN UPPER CASE  !!!!!!!!!!!!!!!! 
degeneracy_tag_CDS <- function(Sequence,Genetic_code = 5, sens = "F",Codon_ = F) {
  library(seqinr)
  library(stringr)
  
  print("Genetic codes are numbered 1-26")
  print("Look at seqinr::translate help to see all the possibilities")
  print("Default is set at Genetic_code = 5, corresponding to the Invertebrate Mitochondrial code")
  print("sens is set to forward as default")
  
  # If given a codon
  if (Codon_ != F) { 
    if (sum(nchar(Sequence)) != 3) { stop("The function is set in codon mode but the given sequence is not 3 characters long") } else {
      if (Codon_ %in% c(1,2,3) == F) { stop("The codon parameter needs to be set at a value between 1 and 3, specifying the position of the nucleotide of interest") } else {
        # If the codon is in one string it needs to be separated into 3 distinc characters
        if (length(Sequence) == 1) { Sequence <- str_split_1(Sequence,pattern="") } else {
          if((length(Sequence) %in% c(1,3)) == F) { stop("The Sequence has to be as un string of 3 characters or a vector of 3 characters") } }
        if (length(Sequence) == 3) {
          
          AA <- seqinr::translate(Sequence)
          nucleotide <- Sequence[Codon_]
          codon <- Sequence
          
          nb_syn <- 0
          for (mut in c("A","T","C","G")){

            if (nucleotide != mut) {
              
              codon_mut <- codon
              codon_mut[Codon_] <- mut
              codon_tr = seqinr::translate(codon,numcode = Genetic_code, sens = sens)
              codon_mut_tr = seqinr::translate(codon_mut,numcode = Genetic_code, sens = sens)
              
              # Count how many times it finds a synonymous mutation
              if (codon_tr == codon_mut_tr) { nb_syn <- nb_syn+1 }
            }
          }
          
          if (nb_syn == 0) { assign("Degeneracy","Zerofold",envir = globalenv()) }
          if (nb_syn == 1) { assign("Degeneracy","Twofold",envir = globalenv()) }
          if (nb_syn == 2) { assign("Degeneracy","Threefold",envir = globalenv()) }
          if (nb_syn == 3) { assign("Degeneracy","Fourfold",envir = globalenv()) }
          
          print("Result in found in the Degeneracy object")
          }
      }
    }
    }
  
  ## If the seq is not said to be a codon
  if (Codon_ == F) { 
    # If the seq is not a list
    if (is.list(Sequence) == F) { 
      Seq_table <- assign("Seq_table", data.frame(Seq = Sequence,Codon_Pos = rep(1:3), Degeneracy = 0), envir = globalenv())
      colnames(Seq_table) = c("Seq","Codon_Pos","Degeneracy") }
    
    # If the list only has one seq, extract it into a table
    if (is.list(Sequence) && length(Sequence) == 1) {
      assign("Seq_table", data.frame(Seq = Sequence, Codon_Pos = rep(1:3), Degeneracy = 0), envir = globalenv())
      colnames(Seq_table) = c("Seq","Codon_Pos","Degeneracy") }
    
    # If the list has multiple sequences, build a new list with all the resulting tables
    if (is.list(Sequence) && length(Sequence) > 1) {
      
      # get list of genes (it takes the seq name)
      gene.list = names(Sequence)
      
      # create receiving list 
      assign("Seq_table_list",Sequence)
      names(Seq_table_list) <- gene.list
      
      # loop to extract seq from list and do the degeneracy test
      for (n in 1:length(Sequence)) {
        
        gene.name = names(Sequence[n])
        Seq_table <- data.frame(Seq = Sequence[[n]],Codon_Pos = rep(1:3), Degeneracy = 0)
        
        ## To tagg the position for their degeneracy we have to check for each position, how many mutations are synonymous
        ## 0 possible synonymous mutation = Zero-fold degeneracy
        ## 1 possible synonymous mutation = Two-fold degeneracy
        ## 2 possible synonymous mutation = Three-fold degeneracy
        ## 3 possible synonymous mutation = Four-fold degeneracy
        for (pos in 1:nrow(Seq_table)) {
          
          ## Get the codon the position is in
          cod_pos <- Seq_table[pos,"Codon_Pos"]
          if (cod_pos == 1) {codon <- c(Seq_table[pos,"Seq"], Seq_table[pos+1,"Seq"], Seq_table[pos+2,"Seq"])}
          if (cod_pos == 2) {codon <- c(Seq_table[pos-1,"Seq"], Seq_table[pos,"Seq"], Seq_table[pos+1,"Seq"])}
          if (cod_pos == 3) {codon <- c(Seq_table[pos-2,"Seq"], Seq_table[pos-1,"Seq"], Seq_table[pos,"Seq"])}
          
          ## Look for the effects of the mutations
          nucleotide <- Seq_table[pos,"Seq"]
          nb_syn = 0
          
          for (mut in c("A","T","C","G")){
            
            if (nucleotide != mut) {
              
              codon_mut <- codon
              codon_mut[cod_pos] <- mut
              codon_tr = seqinr::translate(codon,numcode = Genetic_code, sens = sens)
              codon_mut_tr = seqinr::translate(codon_mut,numcode = Genetic_code, sens = sens)
              
              # Count how many times it finds a synonymous mutation
              if (codon_tr == codon_mut_tr) { nb_syn <- nb_syn+1 }
              
            }
          }
          if (nb_syn == 0) { Seq_table[pos,"Degeneracy"] <- "Zerofold" }
          if (nb_syn == 1) { Seq_table[pos,"Degeneracy"] <- "Twofold" }
          if (nb_syn == 2) { Seq_table[pos,"Degeneracy"] <- "Threefold" }
          if (nb_syn == 3) { Seq_table[pos,"Degeneracy"] <- "Fourfold" }
        }
        Seq_table_list[[n]] <- Seq_table
        
      }
      assign("Seq_table_list",Seq_table_list, envir = globalenv())
      View(Seq_table_list)
      print("The resulting list with all the genes as tables with each position tagged for the level of degeneracy can be found under Seq_table_list")
    }    
    
    ## To tagg the position for their degeneracy we have to check for each position, how many mutations are synonymous
    # (This same loop is repeated above for the case when sequences are in a list)
    ## 0 possible synonymous mutation = Zero-fold degeneracy
    ## 1 possible synonymous mutation = Two-fold degeneracy
    ## 2 possible synonymous mutation = Three-fold degeneracy
    ## 3 possible synonymous mutation = Four-fold degeneracy
    if (is.list(Sequence) == F || is.list(Sequence) && length(Sequence) == 1) {
      for (pos in 1:nrow(Seq_table)) {
        
        ## Get the codon the position is in
        cod_pos <- Seq_table[pos,"Codon_Pos"]
        if (cod_pos == 1) {codon <- c(Seq_table[pos,"Seq"], Seq_table[pos+1,"Seq"], Seq_table[pos+2,"Seq"])}
        if (cod_pos == 2) {codon <- c(Seq_table[pos-1,"Seq"], Seq_table[pos,"Seq"], Seq_table[pos+1,"Seq"])}
        if (cod_pos == 3) {codon <- c(Seq_table[pos-2,"Seq"], Seq_table[pos-1,"Seq"], Seq_table[pos,"Seq"])}
        
        ## Look for the effects of the mutations
        nucleotide <- Seq_table[pos,"Seq"]
        nb_syn = 0
        
        for (mut in c("A","T","C","G")){
          
          if (nucleotide != mut) {
            
            codon_mut <- codon
            codon_mut[cod_pos] <- mut
            codon_tr = seqinr::translate(codon,numcode = Genetic_code, sens = sens)
            codon_mut_tr = seqinr::translate(codon_mut,numcode = Genetic_code, sens = sens)
            
            # Count how many times it finds a synonymous mutation
            if (codon_tr == codon_mut_tr) { nb_syn <- nb_syn+1 }
            
          }
        }
        if (nb_syn == 0) { Seq_table[pos,"Degeneracy"] <- "Zerofold" }
        if (nb_syn == 1) { Seq_table[pos,"Degeneracy"] <- "Twofold" }
        if (nb_syn == 2) { Seq_table[pos,"Degeneracy"] <- "Threefold" }
        if (nb_syn == 3) { Seq_table[pos,"Degeneracy"] <- "Fourfold" }
      }
      assign("Seq_table",Seq_table,envir = globalenv())
      View(Seq_table)
      print("The resulting list with all the genes as tables with each position tagged for the level of degeneracy can be found under Seq_table")
      
    }}
  
  
  
}
####################################################################################################################

####################################################################################################################
### Function to compute Variant frequency of SNPs throught Sanger sequencing
# The script compute average peak height to normalise the peak height at the variant's location
# Inputs are
# file_loc = location of the ab1 file (Sanger seq)
# VAF_pos = position where the variant that wanting to be assessed is located
# VAR_1 = Nucleotide of the first variant (A, T, C, G)
# VAR_1 = Nucleotide of the second variant (A, T, C, G)
Sanger_seq_VAF <- function(file_loc, VAF_pos, VAR_1, VAR_2, discard = 50){
  library(sangerseqR)
  ## Read ab1 file
  file <- read.abif(file_loc)
  
  ##### ab1 file processing #####
  {
    # Generate a sangerseq object (performs base calling)
    seq_obj <- sangerseq(file)
    
    # Extract the base calls
    called_bases <- as.character(primarySeq(seq_obj))
    
    # Get the trace matrix (raw intensity values for A, C, G, T at each trace index)
    traces <- traceMatrix(seq_obj)  # rows = positions, cols = A, C, G, T
    colnames(traces) <- c("A","C","G","T")
    
    # Get the positions in the trace matrix that correspond to each base call
    # (This aligns trace indices to base positions)
    peak_positions <- peakPosMatrix(seq_obj)  # same number of rows as base calls
    colnames(peak_positions) <- c("A","C","G","T")
    
    # Initialize a vector to hold best peak height per base
    best_peak_heights <- numeric(length(called_bases))
    
    # Match trace channel to base call
    channels <- colnames(traces)  # should be "A", "C", "G", "T"
    
    # Split calles bases
    called_bases <- unlist(strsplit(as.character(primarySeq(seq_obj)), split = ""))
    
    # Get approximate peak center positions (e.g., every ~10 trace points for each base)
    # These positions come from the ABIF field PLOC (peak locations)
    base_positions <- file@data$PLOC.2
    
    # Window size for peak search around expected position
    window_half_width <- 5
    
    for (i in seq_along(called_bases)) {
      base <- called_bases[i]
      if (!base %in% c("A", "C", "G", "T")) next
      
      trace_col <- which(channels == base)
      center <- base_positions[i]
      
      if (is.na(center)) next
      
      # Define a small window around the expected peak location
      window_range <- (center - window_half_width):(center + window_half_width)
      window_range <- window_range[window_range > 0 & window_range <= nrow(traces)]  # keep valid indices
      
      if (length(window_range) == 0) next
      
      # Get the max intensity in that window for this base
      best_peak_heights[i] <- max(traces[window_range, trace_col], na.rm = TRUE)
    }
    
    ### Get peak height for each nucleotide at each pos
    # Initialize matrix to hold max peak heights for all bases at each position
    all_peak_heights <- matrix(NA, nrow = length(called_bases), ncol = 4)
    colnames(all_peak_heights) <- c("A_height", "C_height", "G_height", "T_height")
    
    for (i in seq_along(called_bases)) {
      center <- base_positions[i]
      if (is.na(center)) next
      
      window_range <- (center - window_half_width):(center + window_half_width)
      window_range <- window_range[window_range > 0 & window_range <= nrow(traces)]
      
      if (length(window_range) == 0) next
      
      for (j in 1:4) {
        base_letter <- channels[j]
        all_peak_heights[i, j] <- max(traces[window_range, j], na.rm = TRUE)
      }
    }
    
    # Combine seq and peak height (best + per nuc)
    Peak_per_Pos_called <- data.frame(
      Base = called_bases,
      Peak_Height = best_peak_heights,
      all_peak_heights
    )
  }
  
  ### Peak_per_Pos_called has best peak height + height per nucl ###
  
  ##### VAF calculations #####
  {
    ## Compute average peak height per nucleotide
    {
      # Discard the pos
      Peak_calc <- Peak_per_Pos_called[discard:nrow(Peak_per_Pos_called),]
      # Compute average
      A_avg <- mean(Peak_calc$Peak_Height[Peak_calc$Base == "A"])
      C_avg <- mean(Peak_calc$Peak_Height[Peak_calc$Base == "C"])
      G_avg <- mean(Peak_calc$Peak_Height[Peak_calc$Base == "G"])
      T_avg <- mean(Peak_calc$Peak_Height[Peak_calc$Base == "T"])
    }
    
    ## Peak VAR
    Peak_VAR_1 <- Peak_per_Pos_called[VAF_pos,paste(VAR_1,"_height",sep="")]
    Peak_VAR_2 <- Peak_per_Pos_called[VAF_pos,paste(VAR_2,"_height",sep="")]
    
    ## Normalise
    N_Peak_VAR_1 <- Peak_VAR_1 / get(paste(VAR_1,"_avg",sep=""))
    N_Peak_VAR_2 <- Peak_VAR_2 / get(paste(VAR_2,"_avg",sep=""))
    
    ## VAF
    VAF_VAR_1 <- N_Peak_VAR_1 / (N_Peak_VAR_1+N_Peak_VAR_2)
    VAF_VAR_2 <- N_Peak_VAR_2 / (N_Peak_VAR_1+N_Peak_VAR_2)
    
    ## Display results
    print("VAR_1 VAF")
    print(VAF_VAR_1)
    print("VAR_2 VAF")
    print(VAF_VAR_2)
  }
}
####################################################################################################################

