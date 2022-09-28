##Helper R script to create a final_masterVCF.txt 
##Sayeh Gorjifard 2022

library(dplyr)
library(tidyr)
library(readr)


table <- read.table("~/GSHackathon/MasterVCF.txt", header=TRUE)
yeast_genome <- read.csv("~/GSHackathon/yeast_genome.csv")
centromeres <- read.csv("~/GSHackathon/centromeres.csv")
gene_protein <- read_delim("~/GSHackathon/gene_protein.tsv", 
                           delim = "\t", escape_double = FALSE, 
                           col_names = FALSE, trim_ws = TRUE) %>%
  dplyr::rename("REGION" = "X1", "GENE" = "X2", "PROTEIN_LENGTH" = "X3", "UNIPROT" = "X4", "DESCRIPTION" = "X5") %>% select(-DESCRIPTION)
#File with GENE information
SGID <- read.csv("~/GSHackathon/gene_protein_SGID.csv") 

foo <- SGID %>% select(REGION,SGDID)

final <- table %>% 
  left_join(yeast_genome %>% select(-POS),by='CHROM') %>% 
  mutate(Chromosome= gsub("chr","",CHROM)) %>% 
  left_join(centromeres,by="CHROM") %>% left_join(gene_protein,by='REGION') %>%
  dplyr::mutate("AA_POS" = stringr::str_extract(PROTEIN, "([0-9])+")) %>% 
  left_join(foo,by="REGION")
  
  
write.table(final,file = "final_MASTERVCF.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

## making final VCF for all experiments
allVCF <- read.csv("all_yEvo_vcf.csv")


finalVCF <- allVCF %>% 
  left_join(yeast_genome %>% select(-POS),by='CHROM') %>% 
  mutate(Chromosome= gsub("chr","",CHROM)) %>% 
  left_join(centromeres,by="CHROM") %>% left_join(gene_protein %>% select(-GENE),by='REGION') %>%
  dplyr::mutate("AA_POS" = stringr::str_extract(PROTEIN, "([0-9])+")) %>% 
  left_join(foo,by="REGION") %>%
  filter(!is.na(instructor))

write.csv(finalVCF,file = "final_allVCF.csv", row.names = FALSE)

