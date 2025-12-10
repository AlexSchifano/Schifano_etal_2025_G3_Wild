setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")

###### Load the script with the functions I made ######
source("~/Documents/Scripts/R_FUNCTIONS.R")

# List of Hawaiian strains within a same cluster (There are 2 NZ/Au strains, ECA2581 & ECA36 as well as 1 Cali QZ1211)
list.hw.Strains <- c("ECA1240","ECA2291","ECA742","ECA1821","ECA1761","ECA1725","ECA1851"
                     ,"ECA1717","ECA347","ECA1192","ECA1251","ECA36","ECA2581","ECA2360"
                     ,"QX1211","ECA2043","ECA740","ECA1843","ECA724","ECA1253","ECA1185")

list.hw.Node <- c("V687","V684","V683","V681","V678","V677","V675","V673","V670",
                  "V685","V667","V665","V651","V649","V654","V664","V657","V663","V662","V661")

list.hw <- c(list.hw.Node,list.hw.Strains)

## Load mutation table
MEGA <- read.table("R_txt_files/Mutation_List.txt",header=T,sep="\t")

# Subset the hawaiian
MEGA_HW <- MEGA[MEGA$Strain_Name %in% list.hw,]

# Load the reference fasta (needed to count number of bp)
library(seqinr)
V170 <- read.fasta("../01-DATA/9-VARIANT_CALLING/REF/V170.fa")

# Mut spectrum
spectrum_plot(MEGA_HW)
