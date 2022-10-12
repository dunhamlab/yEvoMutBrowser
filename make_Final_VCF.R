##Helper R script to create a final_masterVCF.txt 
##Sayeh Gorjifard 2022

library(dplyr)
library(tidyr)
library(readr)


table <- read.table("MasterVCF.txt", header=TRUE)
yeast_genome <- read.csv("yeast_genome.csv")
centromeres <- read.csv("centromeres.csv")
gene_protein <- read_delim("gene_protein.tsv", 
                           delim = "\t", escape_double = FALSE, 
                           col_names = FALSE, trim_ws = TRUE) %>%
  dplyr::rename("REGION" = "X1", "GENE" = "X2", "PROTEIN_LENGTH" = "X3", "UNIPROT" = "X4", "DESCRIPTION" = "X5") %>% select(-DESCRIPTION)
#File with GENE information
SGID <- read.csv("gene_protein_SGID.csv") 

SGID <- SGID %>% select(REGION,SGDID)

final <- table %>% 
  #we're adding the chromosome lengths to the table by merge on CHROM column
  left_join(yeast_genome %>% select(-POS),by='CHROM') %>% 
  #make a new column for table and take out the "chr" in the CHROM column such that the new column only has the number 
  mutate(Chromosome= gsub("chr","",CHROM)) %>% 
  #add the centromeme 1 and 2 values to table by joining on CHROM column 
  left_join(centromeres,by="CHROM") %>% left_join(gene_protein,by='REGION') %>%
  #changing the insides of the AA_pos column
  dplyr::mutate("AA_POS" = stringr::str_extract(PROTEIN, "([0-9])+")) %>% 
  #left_join SGID by the region column so I have the ID's of the genes for the website
  left_join(SGID,by="REGION")
  
  
write.table(final,file = "final_MASTERVCF.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)


######### the actual master table we used from rene post-hackathon 
## making final VCF for all experiments
allVCF <- read.csv("all_yEvo_vcf.csv")


finalVCF <- allVCF %>% 
  left_join(yeast_genome %>% select(-POS),by='CHROM') %>% 
  mutate(Chromosome= gsub("chr","",CHROM)) %>% 
  left_join(centromeres,by="CHROM") %>% left_join(gene_protein %>% select(-GENE),by='REGION') %>%
  dplyr::mutate("AA_POS" = stringr::str_extract(PROTEIN, "([0-9])+")) %>% 
  left_join(SGID,by="REGION") %>%
  filter(!is.na(instructor))

write.csv(finalVCF,file = "final_allVCF.csv", row.names = FALSE)

