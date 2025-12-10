setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")

###### Making variant freq violin plots #####
{
  ## I have to use the Var table which has the variant info
  # Reading the variant table (do not need the fixed differences)
  Var = read.table("R_txt_files/Variants_List.txt",header = T, sep="\t")
  library(ggplot2)
  library(gcookbook)
  
  ##### Plot for all variants together ######
  {
    print(ggplot(Var, aes(x=0, y=Variant_Freq)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
      theme_minimal() + ylim(-0.05,1.05) +
      labs(y = "Variant Frequency") + ggtitle("No filter All mutations") +
      theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank()))
  }
  
  ##### Plot for all variants colored by gene type #####
  {
    print(ggplot(Var, aes(x=0, y=Variant_Freq,fill=Gene_Type)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
      theme_minimal() + ylim(-0.05,1.05) + ggtitle("No filter Coloured by genetype") +
      labs(y = "Variant Frequency") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank()) )
  }
  
  ##### Plot with removing fixed variants (fixed in the isotypes but not in tree ASR) #####
  {
    ### Setting the threshold at 1 ###
    {
      # Removing the fixed variants
      Var_nofixed1 <- Var[Var$Variant_Freq < 1,]
      
      # Plotting all together
      print(ggplot(Var_nofixed1, aes(x=0, y=Variant_Freq)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
        theme_minimal() + ylim(-0.05,1.05) + ggtitle("VAF < 1 All mutations") +
        labs(y = "Variant Frequency") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank()))
      
      # Plotting with coloring per gene type
      print(ggplot(Var_nofixed1, aes(x=0, y=Variant_Freq,fill=Gene_Type)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
        theme_minimal() + ylim(-0.05,1.05) + ggtitle("VAF < 1 Coloured by genetype") +
        labs(y = "Variant Frequency") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank()) )
    }
    
    ### Setting the threshold at .97 ###
    {
      # Removing the fixed variants
      Var_nofixed97 <- Var[Var$Variant_Freq < .97,]
      
      # Plotting all together
      print(ggplot(Var_nofixed97, aes(x=0, y=Variant_Freq)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
        theme_minimal() + ylim(-0.05,1.05) + ggtitle("VAF < .97 All mutations") +
        labs(y = "Variant Frequency") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank()))
      
      # Plotting with coloring per gene type
      print(ggplot(Var_nofixed97, aes(x=0, y=Variant_Freq,fill=Gene_Type)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
        theme_minimal() + ylim(-0.05,1.05) + ggtitle("VAF < .97 coloured by genetye") +
        labs(y = "Variant Frequency") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank()) )
    }
    
    ### Making a composite plot with all variants and only the heteroplasmic ones ( .03 < VAF < .97) ###
    {
      Var_nofixed97 <- Var[Var$Variant_Freq < .97,]
      
      all_mut <- ggplot(Var, aes(x=0, y=Variant_Freq)) + geom_violin(trim = F,scale = "count", bw = 0.01, fill = "gray70") + 
              theme_minimal() + ylim(-0.05,1.05) + 
              geom_violin(data = Var_nofixed97, aes(x=1, y=Variant_Freq),trim = F,scale = "count", bw=0.01, fill = "cadetblue") +
              labs(y = "") +
              theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank())
      
      ggsave("Writing_NEW/5-Figures/Figures/PDF/Figure6/RAW/All_mut_VAF.pdf", plot = all_mut, height = 20, width = 20)
    
    }
    
  }
  
  ##### Plotting NSYN vs SYN #####
  {
    # I have to use CDS muts only
    Var_CDS <- Var[Var$Gene_Type == "CDS",]
    
    # Plotting with coloring per synonymous
    print(ggplot(Var_CDS, aes(x=0, y=Variant_Freq,fill=Synonymous)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
      theme_minimal() + ylim(-0.05,1.05) + ggtitle("No filter CDS mutations, coloured by Effect") +
      labs(y = "Variant Frequency") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank()))
    
    ## Same for nofixed but threshold at 97% (includes seq error)
    Var_CDS_nofixed97 <- Var_CDS[Var_CDS$Variant_Freq < .97,]
    
    # Plotting with coloring per synonymous
    print(ggplot(Var_CDS_nofixed97, aes(x=0, y=Variant_Freq,fill=Synonymous)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
      theme_minimal() + ylim(-0.05,1.05) + ggtitle("VAF < .97 CDS mutations, coloured by Effect") +
      labs(y = "Variant Frequency") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank())) 
    
    ### Making a composite with both thresholds ###
    {
      library(ggpattern)
      
      CDS_mut <- ggplot() + 
        geom_violin_pattern(data = Var_CDS, aes( x = 0, y = Variant_Freq, group = Synonymous, pattern = Synonymous), 
                            trim = F,scale = "count", bw = 0.01,fill="gray75", pattern_fill = "black",
                            pattern_density = 0.1, pattern_spacing = 0.02) + 
              theme_minimal() + ylim(-0.05,1.05) +
              geom_violin_pattern(data = Var_CDS_nofixed97, aes(x = 1, y = Variant_Freq, group = Synonymous, pattern = Synonymous),
                          trim = F,scale = "count", bw = 0.01,fill="cadetblue", pattern_fill = "black",
                          pattern_density = 0.1, pattern_spacing = 0.02) + 
              labs(y = "") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank())
      
      ggsave(filename = "Writing_NEW/5-Figures/Figures/PDF/Figure6/RAW/CDS_mut_VAF.pdf",plot=CDS_mut,height = 20,width = 20)
    }
      }

  ##### Plotting INDELs #####
  {
    Var_INDEL <- Var[Var$Mut_Type == "INDEL",]
    
    # Plotting with coloring per genetype
    VAF_Indels <- ggplot(Var_INDEL, aes(x=0, y=Variant_Freq,fill=Gene_Type)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
      theme_minimal() + ylim(-0.05,1.05) + ggtitle("No filter INDELs, coloured by genetype") + scale_fill_manual(values = c("darkred","grey90")) +
      labs(y = "Variant Frequency") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank())
    
    ggsave(VAF_Indels,width = 20,height = 20,filename = "Writing_NEW/5-Figures/Figures/PDF/Figure6/RAW/VAF_Indels.pdf")
    
    ##### Remove fixed VAF < 1 #####
    {
      Var_INDEL_nofixed1 <- Var_INDEL[Var_INDEL$Variant_Freq < 1,]
      
      # Plotting with coloring per genetype
      print(ggplot(Var_INDEL_nofixed1, aes(x=0, y=Variant_Freq,fill=Gene_Type)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
        theme_minimal() + ylim(-0.05,1.05) + ggtitle("VAF < 1 INDELs, coloured by genetype") +
        labs(y = "Variant Frequency") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank())) 
    }
    
    ##### Remove fixed VAF < .97 #####
    {
      Var_INDEL_nofixed97 <- Var_INDEL[Var_INDEL$Variant_Freq < .97,]
      
      # Plotting with coloring per genetype
      print(ggplot(Var_INDEL_nofixed97, aes(x=0, y=Variant_Freq,fill=Gene_Type)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
        theme_minimal() + ylim(-0.05,1.05) + ggtitle("VAF < .97 INDELs, coloured by genetype") +
        labs(y = "Variant Frequency") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank())) 
    }

  }
}

##### Looking at proportions of each variant frequency classes #####
{
  # Fixed variants
  Fixed <- Var[Var$Variant_Freq == 1,]
  Fixed_prop <- nrow(Fixed)/nrow(Var)
  
  # All variants above heteroplasmy treshold (kind of fixed)
  HighFreq <- Var[Var$Variant_Freq >.97,]
  HighFreq_prop <- nrow(HighFreq)/nrow(Var)
  
  # Looking at only heteroplasmic variants (.03 < VAF < 0.97)
  # There are no variants below .07
  Var_nofixed97 <- Var[Var$Variant_Freq < .97,]
  hist(Var_CDS_nofixed97$Variant_Freq)
  nrow(Var_nofixed97[Var_nofixed97$Variant_Freq > .79,]) / nrow(Var_nofixed97)
  nrow(Var_nofixed97[Var_nofixed97$Variant_Freq < .79,]) / nrow(Var_nofixed97)
  
  # How many variants between .4 and .9 VAF
  Var_VAF4 <- Var[Var$Variant_Freq > .4,]
  Var_VAF49 <- Var_VAF4[Var_VAF4$Variant_Freq < .9,]
  nrow(Var_VAF49) / nrow(Var_CDS_nofixed97)
}

##### Konrad Variant data #####
{
  Konrad <- read.table("Previous_results/Konrad/Table1_Konrad2017_MALines_raw.txt", fill=T,header=T, sep="\t")
  Konrad$Frequency <- as.numeric(Konrad$Frequency)
  
  # Variant freq plot
  Konrad_VAF <- ggplot(Konrad, aes(x=0, y=Frequency)) + geom_violin(trim = F,scale = "count", bw = 0.02) + 
          theme_minimal() + ylim(-0.05,1.05) +
          labs(y = "Variant Frequency") + 
          theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank())
  
  ggsave(Konrad_VAF,width = 20,height = 20,filename = "Writing_NEW/5-Figures/Figures_pptx/Figure6/VAF_Konrad.pdf")

  }


##### Removing multiple heteroplasmies
{
  Var_nofixed97 <- Var[Var$Variant_Freq < .97,]
  tab_Hetero <- data.frame(table(Var_nofixed97$Strain_Name))
  tab_Hetero_2 <- tab_Hetero[tab_Hetero$Freq == 1,]
  
  Var_nofixed97_RM <- Var_nofixed97[(Var_nofixed97$Strain_Name %in% tab_Hetero_2$Var1),]
  
  print(ggplot(Var_nofixed97_RM, aes(x=0, y=Variant_Freq)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
          theme_minimal() + ylim(-0.05,1.05) + ggtitle("VAF < .97 removing several heteroplasmy lines") +
          labs(y = "Variant Frequency") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank())) 
  
  print(ggplot(Var_nofixed97_RM, aes(x=0, y=Variant_Freq,fill=Synonymous)) + geom_violin(trim = F,scale = "count", bw = 0.01) + 
          theme_minimal() + ylim(-0.05,1.05) + ggtitle("VAF < .97 removing several heteroplasmy lines") +
          labs(y = "Variant Frequency") + theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank())) 
}
