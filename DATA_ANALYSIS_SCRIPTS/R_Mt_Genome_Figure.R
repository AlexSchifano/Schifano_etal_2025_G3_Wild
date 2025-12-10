setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/Scripts/")

# install.packages("circlize")
library(circlize)

#### Create a Data Frame for the *C. elegans* mtDNA Genes ####
{
  genes <- data.frame(
    gene = c("tRNA-P", "NC1", "tRNA-V", "nd-6","NC2", "nd-4l","NC3", "tRNA-W", "tRNA-E", "rRNA-12s", "tRNA-S1", "tRNA-N","NC4", "tRNA-Y", 
             "nd-1", "atp-6","NC5", "tRNA-K", "tRNA-L1", "tRNA-S2", "nd-2","NC6", "tRNA-I", "tRNA-R", "tRNA-Q","NC7", "tRNA-F", 
             "ctb-1","NC8", "tRNA-L2", "cox-3","NC9", "tRNA-T", "nd-4","NC10", "cox-1", "tRNA-C", "tRNA-M", "tRNA-D", "tRNA-G", 
             "cox-2","NC11", "tRNA-H", "rRNA-16s", "nd-3", "nd-5","tRNA-A","D-loop"),
    
    start = c(1, 56, 58, 113, 548, 549, 783, 785, 842, 898, 1595, 1648, 1704, 1707, 1763, 2634, 3234, 3244, 3307, 3362, 3418,
              4267, 4269, 4330, 4385, 4440, 4447, 4504, 5617, 5621, 5678, 6446, 6450, 6506, 7736, 7845, 9422, 9478, 9538, 9593, 
              9649, 10345, 10348, 10403, 11356, 11691, 13275, 13329),
    
    end = c(55, 57, 112, 547, 548, 782, 784, 841, 897, 1594, 1647, 1703, 1706, 1762, 2638, 3233, 3243, 3306, 3361, 3417, 4266,
            4268, 4329, 4384, 4439, 4446, 4503, 5616, 5620, 5677, 6445, 6449, 6505, 7735, 7844, 9422, 9477, 9537, 9592, 
            9648, 10344, 10347, 10402, 11355, 11691, 13274, 13328, 13794),
    
    color = "gray50"
  )
  
  # Making type groupings
  tRNA_list <- grep("tRNA-",x = genes$gene)
  rRNA_list <- grep("rRNA-",x = genes$gene)
  ETC1 <- grep("nd",x = genes$gene)
  ETC3 <- grep("ctb-",x = genes$gene)
  ETC4 <- grep("cox-",x = genes$gene)
  ETC5 <- grep("atp-",x = genes$gene)
  
  # Adding a color for each group, gray50 will be the color of the rest = non coding + D-loop
  genes[tRNA_list,"color"] <- "darkred"
  genes[rRNA_list,"color"] <- "tan1"
  genes[ETC1,"color"] <- "steelblue3"
  genes[ETC3,"color"] <- "paleturquoise1"
  genes[ETC4,"color"] <- "purple"
  genes[ETC5,"color"] <- "yellowgreen"
}

library(circlize)

genes <- genes[grep(x =genes$gene, pattern = "NC", invert = T),]
genes <- genes[grep(x =genes$gene, pattern = "tRNA", invert = T),]
xlim <- data.frame(as.numeric(genes$start), as.numeric(genes$end))

sectors <- factor(genes$gene, levels = genes$gene)

circos.initialize(sectors, xlim = xlim)
circos.track(genes$gene, ylim = c(0,1))
circos.info(plot=T,sector.index = F)
