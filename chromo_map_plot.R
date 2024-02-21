# Load the required libraries
library(ggplot2)
library(plotly)
library(dplyr)
library(forcats)

# Load data from a CSV file
chromosometest <- read.csv(file.path(getwd(), "/chromosome_info.csv"))
genetest <- read.csv(file.path(getwd(), "/gene_info.csv"))
mutations_test <- read.csv(file.path(getwd(),"/all_yEvo_vcf.csv"))

#test 150 genes
#genetest <- read.csv("C:/Users/virgi/OneDrive/UW/Other/test-yEvo/gene_info_test.csv")

filtered_mutations <-filter(mutations_test, condition == "caffeine")
# head(mutations_filtered)
# Assuming your dataframes are mutations_filtered and gene_test

# Merge the dataframes based on the "REGION" column
filtered_genes_data <- left_join(filtered_mutations, genetest, by = "REGION")
filtered_genes_data <- filtered_genes_data[complete.cases(filtered_genes_data), ]

# Define the desired order of categories
desired_order <- c('chrM', 'chrXVI', 'chrXV', 'chrXIV', 'chrXIII', 'chrXII', 'chrXI', 'chrX', 'chrIX', 'chrVIII', 'chrVII', 'chrVI', 'chrV', 'chrIV', 'chrIII', 'chrII', 'chrI')


# Convert category to a factor with the desired order
chromosometest$CHROM <- factor(chromosometest$CHROM, levels = desired_order)

# Sample data for chromosomes
chromosomes <- data.frame(
  chromosome = chromosometest$CHROM,
  length = chromosometest$length
)

# Sample data for genes
genes <- data.frame(
  chromosome = genetest$CHROM,
  start = genetest$START,
  end = genetest$STOP,
  geneName = genetest$GENE 
)


filtered_genes <- data.frame(
  chromosome = filtered_genes_data$CHROM.y,
  start = filtered_genes_data$START,
  end = filtered_genes_data$STOP,
  geneName = filtered_genes_data$GENE.x
)

library(dplyr)

# Create an empty dataframe to store the final results
final_gene <- data.frame()

# Iterate through unique genes
unique_genes <- unique(filtered_genes_data$GENE.x)
for (gene in unique_genes) {
  # Filter data for the current gene
  gene_data <- filtered_genes_data %>% filter(GENE.x == gene)
  
  # Count the number of repetitions
  num_repeats <- nrow(gene_data)
  
  # Get the associated columns
  gene_info <- gene_data[1, c("CHROM.x", "START", "STOP", "GENE.x")]
  
  # Add the number of repeats as a new column
  gene_info$repeats <- num_repeats
  
  # Add this gene to the final dataframe
  final_gene <- rbind(final_gene, gene_info)
}

# Print the final dataframe
print(final_gene)





# Define a mapping from chromosome names to numbers
chromosome_mapping <- c(
  'chrM' = 1, 
  'chrXVI' = 2, 
  'chrXV' = 3, 
  'chrXIV' = 4, 
  'chrXIII' = 5, 
  'chrXII' = 6, 
  'chrXI' = 7, 
  'chrX' = 8, 
  'chrIX' = 9, 
  'chrVIII' = 10, 
  'chrVII' = 11, 
  'chrVI' = 12, 
  'chrV' = 13, 
  'chrIV' = 14, 
  'chrIII' = 15, 
  'chrII' = 16, 
  'chrI' = 17
)

# Add a new column to the dataframe with mapped chromosome numbers
filtered_genes$chromosome_as_num <- chromosome_mapping[filtered_genes$chromosome]
chromosomes$chromosome_as_num <- chromosome_mapping[chromosomes$chromosome]
final_gene$chromosome_as_num <- chromosome_mapping[final_gene$CHROM.x]

# -------------------------------------------------------------------------



# Display the updated dataframe
print(filtered_genes)
print(chromosomes)

# Plotting
p <- ggplot() +
  geom_bar(data = chromosomes, aes(x = length, y = chromosome), stat = 'identity', fill = 'lightblue', width = 0.5) +

  # geom_rect(data = filtered_genes, aes(ymin = chromosome_as_num - 0.4,
  #                                      ymax = chromosome_as_num + 0.4,
  #                                      xmin = start, 
  #                                      xmax = end,
  #                                      text = geneName), 
  geom_rect(data = final_gene, aes(ymin = chromosome_as_num - 0.4,
                                       ymax = chromosome_as_num + 0.4,
                                       xmin = START, 
                                       xmax = STOP,
                                       text = paste("Gene Name: ",GENE.x),
                                       fill = repeats), 
          
           , alpha = 0.5) +
  scale_fill_gradient(low = "navy", high = "red",limits = c(0, 25)) +
  #geom_rect(data = filtered_genes, aes(ymin = 1, ymax = 20, xmin = 1, xmax = 20), fill = "red") + 
  labs(title = 'Chromosomes with Genes Overlay',
       y = 'Chromosome',
       x = 'Length')
# Print the plot
print(p)

# Convert ggplot2 plot to plotly
ggplotly(p)
print("done")
