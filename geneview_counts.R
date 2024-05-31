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
  
# Iterate through unique genes
count_proteins_same <- cur_gene %>%
  group_by(POS) %>%
  summarize(
    GENE = first(GENE),
    PROTEIN = first(PROTEIN),
    ANNOTATION = first(ANNOTATION),
    Counts = n(),
    Letter1 = substr(PROTEIN, 1, 1),  # Extract the first character
    Numbers = as.numeric(str_extract(PROTEIN, "[0-9]+")),  # Extract the numbers
    Letter2 = substr(PROTEIN, nchar(PROTEIN), nchar(PROTEIN))
  ) %>%
  ungroup()

count_proteins_same2 <- count_proteins_same %>%
  group_by(Numbers) %>%
  summarize(
    GENE = first(GENE),
    PROTEIN = paste((PROTEIN), collapse = ", "),
    ANNOTATION = first(ANNOTATION),
    list_counts = paste((paste(Counts)), collapse = ", "),
    Counts_tot = sum(Counts),
  ) %>%
  ungroup()

count_proteins_same2$PROTEIN <- sapply(count_proteins_same2$PROTEIN, FUN = strsplit, split = ",", simplify = TRUE)
count_proteins_same2$list_counts <- sapply(count_proteins_same2$list_counts, FUN = strsplit, split = ",", simplify = TRUE)

# test <- rbind((count_proteins_same2$PROTEIN), (count_proteins_same2$list_counts))
for (i in 1:nrow(count_proteins_same2)) {
  # Access row elements using row index (i)
  # Or access specific columns
  current_x <- list(count_proteins_same2$list_counts[i])

  current_y <- list(count_proteins_same2$PROTEIN[i])

  my_list <- list()
  print("current_x")
  print(current_x)
  for (j in 0:length(current_x)) {
    print("current_x[j]")
    test <- current_x[j]
    print(test)
    string_test <- paste(test , ": " ,(current_y[j]))
    my_list <- c(my_list,string_test )
    # print('string test')
    # print(string_test)
  }
  # print('full list')
  # print(my_list)
  
  # count_proteins_same2$help[i] <- rbind(list(current_y), list(current_x))
}

print('TESTETSTET')
current_protein2 <- count_proteins_same2$PROTEIN[3]
# test_protein <- strsplit(current_protein2, ",")[[1]]
print(current_protein2)
print(typeof(current_protein2))
print(length(current_protein2))
print(list(current_protein2))
print(length(list(current_protein2)))
split_string <- strsplit(current_protein2[[1]], ",")
print(split_string)
print(length(split_string))
# Access elements (assuming no spaces within elements)
element1 <- split_string[1]
element2 <- split_string[2]



combined_strings <- character(nrow(count_proteins_same2))  # Pre-allocate character vector
for (i in 1:nrow(count_proteins_same2)) {
  current_protein <- (count_proteins_same2$PROTEIN[i])
  current_counts <- count_proteins_same2$list_counts[i]

  split_string <- strsplit(current_protein[[1]], ",")
  split_counts <- strsplit(current_counts[[1]], ",")
  
  if (length(split_string) > 1) {
    cur <- list()
    for (j in 1:length(split_string)) {
      current_protein = split_string[j]
      current_protein <- str_trim(current_protein)
      current_counts = split_counts[j]
      Letter1 <- substr(current_protein, 1, 1)  # Extract the first character
      Numbers <- as.numeric(str_extract(current_protein, "[0-9]+"))  # Extract the numbers
      Letter2 <- substr(current_protein, nchar(current_protein), nchar(current_protein))
      print(current_protein)
      print(Letter1)
      cur <- c(cur, paste("Count ", Letter1, '->', Letter2, ": ", current_counts, "\n"))
    }
    combined_strings[i] <- sapply(cur, function(x) paste(x)) %>% paste(collapse = ", ")
    }
  else {
  Letter1 <- substr(current_protein, 1, 1)  # Extract the first character
  Numbers <- as.numeric(str_extract(current_protein, "[0-9]+"))  # Extract the numbers
  Letter2 <- substr(current_protein, nchar(current_protein), nchar(current_protein))
  combined_strings[i] <- paste("Count ", Letter1, '->', Letter2, ": ", current_counts)
  }
}

# Add the combined strings as a new column
count_proteins_same2$combined <- combined_strings







# Define a function to separate and count mutations
separate_counts <- function(x) {
  mutations <- strsplit(x, ", ")[[1]]  # Split protein string by comma
  counts <- strsplit(x, ", ")[[2]] %>% as.numeric()  # Split counts and convert to numeric
  paste(sapply(mutations, paste, collapse = ": "), collapse = ", ")  # Combine mutations and counts
}







