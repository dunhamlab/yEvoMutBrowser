library(dplyr)
library(stringr)

mut_backend <- read.csv("all_yEvo_vcf.csv")
genes_info <- read.csv("gene_info.csv")
chromosome_mapping <- c(
  chrM = 1, chrXVI = 2, chrXV = 3, chrXIV = 4, chrXIII = 5, chrXII = 6,
  chrXI = 7, chrX = 8, chrIX = 9, chrVIII = 10, chrVII = 11, chrVI = 12,
  chrV = 13, chrIV = 14, chrIII = 15, chrII = 16, chrI = 17
)
# Create an empty dataframe to store the final results

gene_name <- "PDR1"

mutation_data_value <- mut_backend

# Merge the data frames based on the “REGION” column
common_cols <- intersect(colnames(mutation_data_value), colnames(genes_info))
mutation_data_value <- merge(mutation_data_value, genes_info, by = common_cols)

mutation_data_value <- mutation_data_value[order(mutation_data_value$GENE),]

# Using subset() function
cur_gene <- subset(mutation_data_value, GENE == gene_name)

# Iterate through unique genes
count_proteins <- cur_gene %>%
    group_by(POS) %>%
  summarize(
      GENE = first(GENE),
      PROTEIN = first(PROTEIN),
      ANNOTATION = first(ANNOTATION),
      Counts = n(),
    ) %>%
  ungroup()
   

# Sample dataframe
df <- data.frame(PROTEIN = c("R324F", "T542R", "A123B"))

# Extracting components using mutate and str_extract
df <- df %>%
  mutate(Letter1 = substr(PROTEIN, 1, 1),  # Extract the first character
         Numbers = as.numeric(str_extract(PROTEIN, "[0-9]+")),  # Extract the numbers
         Letter2 = substr(PROTEIN, nchar(PROTEIN), nchar(PROTEIN)))  # Extract the last character

# Display the updated dataframe
print(df)
