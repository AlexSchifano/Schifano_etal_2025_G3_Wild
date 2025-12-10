setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")
  
##### McDonald–Kreitman test #####
{
  ## Test for selection by comparing variation between species (substitutions) and within species (polymorphisms) 
  ## I have to look at the protein coding sites only
  library(ape)
  library(Biostrings)
  library(seqinr)
  
  #### Within species fixed polymorphisms #### 
  {
    ASR = read.table("../01-DATA/8-COUNT_MUTATIONS/Count_Mut_results_REF_POS.txt", header = T, sep="\t")
    ASR_CDS <- ASR[ASR$Gene_Type == "CDS",]
    
    ## Get the within species fixed mut counts
    Within_Syn = nrow(ASR_CDS[ASR_CDS$Synonymous == "Synonymous",])
    Within_NSyn = nrow(ASR_CDS[ASR_CDS$Synonymous == "NON_Synonymous",])
  }
  
  #### Between species mut count and annotation ####
  {
    # Read a fasta from Brenneri that has the same trimming as the wild isolates + is only the CDS
    Brenneri <- read.fasta("MKT/Brenneri_CDS_only_string_noATP6stop.fasta")
    
    # Read a fasta of the wild isolates with only the CDS
    # The wild_msa has been made with consensus sequences so it does not include variants information
    wild_CDS_string <- read.fasta("MKT/Wild_MSA_CDS_string.fasta")
    
    ## Counting the mutations between Brenneri and every wild isolate
    
    # Getting the seq in a data frame
    Brenneri <- data.frame(seq(1:length(Brenneri[[1]])),Brenneri[[1]])
    colnames(Brenneri) <- c("pos","seq")
    
    ## Loop for mut detection
    # Loops for each isolate one by one and looks at every mutation that differs from Brenneri
    # Stores every result in a file so delete the file before running the loop
    # The loop prints how many isolates have been processed
    for (strain in names(wild_CDS_string)) {
      
      # Preparing the table of mutations
      mut1=data.frame(matrix(ncol=5,nrow=1))
      colnames(mut1) = c("POS","Mut_from","->","Mut_to","Strain")
      
      # extract elegans wild isolate
      ele <- data.frame(seq(1:length(wild_CDS_string[[strain]])),wild_CDS_string[[strain]])
      colnames(ele) <- c("pos","seq")
      
      # Detect mutations
      x=1
      for (pos in 1:nrow(Brenneri)) {
        if (Brenneri[pos,"seq"] != ele[pos,"seq"]) {
          mut1[x,"POS"] <- pos
          mut1[x,"Mut_from"] <- Brenneri[pos,"seq"]
          mut1[x,"->"] <- "->"
          mut1[x,"Mut_to"] <- ele[pos,"seq"]
          mut1[x,"Strain"] <- strain
          x=x+1
        }}
      print(strain)
      print(which(names(wild_CDS_string) == strain))
      write.table(mut1, file = "R_txt_files/MKT_Count_mut_between.txt", sep = "\t", quote = F, row.names = F, col.names = F, append = T)}
    
    # Read the mut counts file
    btw = read.table("R_txt_files/MKT_Count_mut_between.txt",sep="\t",header = T)
    colnames(btw) <- c("POS","Mut_from","->","Mut_to","Strain")
    
    # Get nb of line per pos found
    pos_count <- data.frame(table(btw$POS))
    colnames(pos_count) <- c("pos","count")
    write.table(pos_count,file = "R_txt_files/MKT_pos_count_btw.txt",quote = F,sep = "\t",row.names = F)
    pos_count <- read.table("R_txt_files/MKT_pos_count_btw.txt",header = T)
  
    # subseting only the positions that have been found in all wild mitotypes
    wild_lines <- length(wild_CDS_string)
    pos_count_shared <- pos_count[pos_count$count == wild_lines,]
    
    ## Making sure the mutations are exactly the same ones, some fixed differences can be polymorphisms in elegans
    # Subsetting the list of mutations to get their details
    btw_mut <- btw[,1:4]
    # Getting the ones at pos shared by all wilds (compared to Brenneri)
    btw_shared <- btw_mut[btw_mut$POS %in% pos_count_shared$pos,]
    # Count the number of variants for each of those positions
    btw_shared_unique <- unique(btw_shared)
    btw_shared_unique_count <- data.frame(table(btw_shared_unique$POS))
    nrow(btw_shared_unique_count[btw_shared_unique_count$Freq == 1,])
    nrow(btw_shared_unique_count[btw_shared_unique_count$Freq == 2,])
    nrow(btw_shared_unique_count[btw_shared_unique_count$Freq == 3,])
    # ==> There are 846 positions that are fixed differences between elegans and brenneri
    # ==> There are 188 fixed between and 2 variants within elegans
    # ==> There are 4 fixed between and 3 variants within elegans
    
    ## Annotating the variants found
    # Preparing the table
    Annotated <- btw_shared_unique
    Annotated$Codon_pos <- ""
    Annotated$Codon_Ele <- ""
    Annotated$Codon_Bre <- ""
    Annotated$AA_Ele <- ""
    Annotated$AA_Bre <- ""
    Annotated$Synonymous <- ""
    
    # Getting the sequences used for translation
    Brenneri <- read.fasta("MKT/Brenneri_CDS_only_string_noATP6stop.fasta")
    Brenneri <- data.frame(seq(1:length(Brenneri[[1]])),Brenneri[[1]])
    colnames(Brenneri) <- c("pos","seq")
    N2 <- read.fasta("MKT/N2_CDS_ref_noATP6stop.fasta")
    N2 <- data.frame(seq(1:length(N2[[1]])),N2[[1]])
    colnames(N2) <- c("pos","seq")
    
    # Annotating 
    for (pos in Annotated$POS) {
      if (pos/3 == ceiling((pos)/3)) {
        # Codon pos
        cod_pos <- 3
        # Codon in Brenneri
        codon_bren <- paste(Brenneri[(pos-2):pos,"seq"],collapse="")
        # AA in Brenneri
        AA_bren <- seqinr::translate(paste(Brenneri[(pos-2):pos,"seq"]),numcode = 5)
        # Codon in N2
        codon_elegans <- paste(N2[(pos-2):(pos),"seq"],collapse="")
        # AA in Elegans
        AA_elegans <- seqinr::translate(paste(N2[(pos-2):(pos),"seq"]),numcode = 5) } 
      if((pos+1)/3 == ceiling((pos)/3)) {
        cod_pos <- 2    
        codon_bren <- paste(Brenneri[(pos-1):(pos+1),"seq"],collapse="")
        AA_bren <- seqinr::translate(paste(Brenneri[(pos-1):(pos+1),"seq"]),numcode = 5)
        codon_elegans <- paste(N2[(pos-1):(pos+1),"seq"],collapse="")
        AA_elegans <- seqinr::translate(paste(N2[(pos-1):(pos+1),"seq"]),numcode = 5) }
      if((pos+2)/3 == ceiling((pos)/3)) {
        cod_pos <- 1
        codon_bren <- paste(Brenneri[pos:(pos+2),"seq"],collapse="")
        AA_bren <- seqinr::translate(paste(Brenneri[pos:(pos+2),"seq"]),numcode = 5)
        codon_elegans <- paste(N2[pos:(pos+2),"seq"],collapse="")
        AA_elegans <- seqinr::translate(paste(N2[pos:(pos+2),"seq"]),numcode = 5) }
      # Adding info to final table
      Annotated$Codon_pos[Annotated$POS == pos] <- cod_pos
      Annotated$Codon_Ele[Annotated$POS == pos] <- codon_elegans
      Annotated$Codon_Bre[Annotated$POS == pos] <- codon_bren
      Annotated$AA_Ele[Annotated$POS == pos] <- AA_elegans
      Annotated$AA_Bre[Annotated$POS == pos] <- AA_bren
      if(AA_elegans == AA_bren) {Annotated$Synonymous[Annotated$POS == pos] <- "Synonymous"} else 
      {Annotated$Synonymous[Annotated$POS == pos] <- "Non_Synonymous"}
    }
    
    ## Add the gene info in the Annotated table
    # Read the GTF (it is made based on the CDS string, the positions are not the ones in the official GTF)
    Brenneri_GTF_CDS_String <- read.table("MKT/Brenneri_CDS_String.gtf",col.names = c("CHRM","START","END","Gene_ID"))
    # Loop to add the gene info in the Annotated table
    for (ro in 1:nrow(Annotated)) {
      pos <- Annotated[ro,"POS"]
      for (g in 1:nrow(Brenneri_GTF_CDS_String)) {
        gene <- Brenneri_GTF_CDS_String[g,"Gene_ID"]
        range <- seq(from = Brenneri_GTF_CDS_String[g,"START"], to = Brenneri_GTF_CDS_String[g,"END"])
        if (pos %in% range) { Annotated[ro,"Gene_ID"] <- gene ; break }
      }
    }
    
    # Writing the annotated file
    write.table(Annotated, file="R_txt_files/MKT_annotated_brenneri_vs_wild_muts.txt",sep = "\t", row.names = F,quote = F)
  }
  # ==> Mutations are in Annotated, "R_txt_files/MKT_annotated_brenneri_vs_wild_muts.txt"
  
  #### Performing the MKT ####
  {
    ## Reading the mutation lists
    ASR = read.table("../01-DATA/8-COUNT_MUTATIONS/Count_Mut_results_REF_POS.txt", header = T, sep="\t")
    ASR_CDS <- ASR[ASR$Gene_Type == "CDS",]
    Annotated = read.table("R_txt_files/MKT_annotated_brenneri_vs_wild_muts.txt", header = T, sep = "\t")
    
    ## Disentangling the positions that are fixed between but polymorphic within
    POS_COUNTS <- data.frame(table(Annotated$POS))
    pos_poly <- POS_COUNTS$Var1[POS_COUNTS$Freq > 1]
    Fixed_poly <- Annotated[Annotated$POS %in% pos_poly,]
    
    ## Fetching the various values
    # DS is number of synonymous muts between species, not counting several times the polymorphic sites
    DS = length(unique(Annotated$POS[Annotated$Synonymous == "Synonymous"]))
    # DN is number if nonsynonymous muts between species
    DN = length(unique(Annotated$POS[Annotated$Synonymous == "Non_Synonymous"]))
    # PS is number of synonymous polymorphic sites (within C. elegans)
    PS = length(unique(ASR_CDS$POS_Full[ASR_CDS$Synonymous == "Synonymous"]))
    # PN is number of nonsynonymous polymorphic sites (within C. elegans)
    PN = length(unique(ASR_CDS$POS_Full[ASR_CDS$Synonymous == "NON_Synonymous"]))
    
    ## Neutrality index
    NI = (PN/PS) / (DN/DS)
    
    ## Test
    table <- data.frame(Poly = c(PN,PS), Fixed = c(DN,DS))
    ftest <- fisher.test(table)
    library(AMR)
    g.test(table)
    chisq.test(table)
    
    #### Doing it gene by gene ####
    {
      # Get gene list to go through
      gene.list <- unique(Annotated$Gene_ID)
      # Prepare the receptacle table
      MKT_Gene <- data.frame(Gene=gene.list,DS="",DN="",PS="",PN="",NI="")
      row.names(MKT_Gene) <- gene.list
      for (ge in gene.list){
        fixed_ge <- Annotated[Annotated$Gene_ID == ge,]
        poly_ge <- ASR_CDS[ASR_CDS$Gene_ID == ge,]
        
        DS_ge <- length(unique(fixed_ge$POS[fixed_ge$Synonymous == "Synonymous"]))
        DN_ge <- length(unique(fixed_ge$POS[fixed_ge$Synonymous == "Non_Synonymous"]))
        PS_ge <- length(unique(poly_ge$POS_Full[poly_ge$Synonymous == "Synonymous"]))
        PN_ge <- length(unique(poly_ge$POS_Full[poly_ge$Synonymous == "NON_Synonymous"]))
        
        NI_ge <- (PN_ge/PS_ge) / (DN_ge/DS_ge)
        MKT_Gene[ge,] <- c(ge,DS_ge,DN_ge,PS_ge,PN_ge,NI_ge) }
      
      MKT_Gene["Total",] <- c("Total",sum(as.numeric(MKT_Gene$DS)),sum(as.numeric(MKT_Gene$DN)),
                              sum(as.numeric(MKT_Gene$PS)),sum(as.numeric(MKT_Gene$PN)),"")
      MKT_Gene["Total","NI"] <- (as.numeric(MKT_Gene["Total","PN"]) / as.numeric(MKT_Gene["Total","PS"])) /
        (as.numeric(MKT_Gene["Total","DN"]) / as.numeric(MKT_Gene["Total","DS"]))
      
      # Writing this table
      write.table(MKT_Gene, file ="R_txt_files/MKT_Per_Gene_Counts.txt", sep = "\t", row.names = F, quote = F )
      
      # Reading the table
      MKT_Gene <- read.table("R_txt_files/MKT_Per_Gene_Counts.txt",header = T)
      row.names(MKT_Gene) <- MKT_Gene$Gene
      
      ## Fisher test on each gene
      for(ge in MKT_Gene$Gene){
        table <- data.frame(Poly = c(MKT_Gene[ge,"PN"],MKT_Gene[ge,"PS"]), Fixed = c(MKT_Gene[ge,"DN"],MKT_Gene[ge,"DS"]))
        ftest <- fisher.test(table)
        chitest <- chisq.test(table,simulate.p.value = T, B = 10000)
        MKT_Gene[ge,"Fsher_PVAL"] <- ftest$p.value
        MKT_Gene[ge,"chisq_10k_PVAL"] <- chitest$p.value
        MKT_Gene[ge,"X-squared"] <- chitest$statistic }
      
      # Writing the final table
      write.table(MKT_Gene, file ="R_txt_files/MKT_Per_Gene_TESTS.txt", sep = "\t", row.names = F, quote = F )
      # also in excel format
      library(openxlsx)
      write.xlsx(MKT_Gene,file = "R_txt_files/MKT_Per_Gene_TESTS.xlsx")
    }
  }
}

  