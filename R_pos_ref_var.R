args = commandArgs()
POS = args[6]
STRAIN = args[7]
ANC = args[8]
POS_REF_CALC = args[9]
TABLE = args[10]

library(stringr)

print(POS)
print(STRAIN)
print(ANC)
print(POS_REF_CALC)
print(TABLE)

## Testing ##
#POS="6739"
#STRAIN="XZ2211"
#ANC="V607"
#POS_REF_CALC="6752"
#TABLE="/Users/alesc399/Documents/C.elegans_sync/Final_files/01-DATA/9-VARIANT_CALLING/ANNOTATION/DATA_Variants_GATK.txt"
##############

# Setting as numeric for matching
POS <- as.numeric(POS)

# Creating the new column if needed
Muts <- read.table(TABLE, sep="\t", header = T)
if (("POS_REF" %in% colnames(Muts)) == F) { Muts$POS_REF <- "" }

# Looking for the matching line and adding the POS_REF
for (ro in 1:nrow(Muts)) {
  if (sum(Muts[ro,c(1,8,9)] == c(POS,STRAIN,ANC)) == 3) { 
    # Add the REF_POS
    Muts[ro,"POS_REF"] <- POS_REF_CALC
    # Add the line to a new file
    LINE <- Muts[ro,]
    write.table(LINE, file = str_replace(TABLE, pattern = ".txt", replacement = "_REF_POS.txt"), quote = F, append = T, row.names = F, col.names = F, sep = "\t")
    }
}

# 