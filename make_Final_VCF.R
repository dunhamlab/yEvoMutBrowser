##Helper R script to create a final_masterVCF.txt 
##Sayeh Gorjifard 2022

library(dplyr)
library(tidyr)

table <- read.table("~/GSHackathon/MasterVCF.txt", header=TRUE)
yeast_genome <- read.csv("/Users/sgorji/Downloads/yeast_genome.csv")
centromeres <- read.csv("~/Documents/centromeres.csv")

final <- table %>% 
  left_join(yeast_genome %>% select(-POS),by='CHROM') %>% 
  mutate(Chromosome= gsub("chr","",CHROM)) %>% 
  left_join(centromeres,by="CHROM")

write.table(final,file = "final_MASTERVCF.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
