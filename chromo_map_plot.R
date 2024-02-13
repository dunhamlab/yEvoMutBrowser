# Load the required libraries
library(ggplot2)
library(plotly)
library(dplyr)

# Load data from a CSV file
chromosometest <- read.csv(file.path(getwd(), "/chromosome_info.csv"))
genetest <- read.csv(file.path(getwd(), "/gene_info.csv"))
mutations_test <- read.csv(file.path(getwd(),"/final_allVCF.csv"))
#test 150 genes
#gene_test <- read.csv("C:/Users/virgi/OneDrive/UW/Other/test-yEvo/gene_info_test.csv")

mutations_filtered <-subset(mutations_test, instructor == "Moscow")
# head(mutations_filtered)


# Assuming your dataframes are mutations_filtered and gene_test

# Merge the dataframes based on the "REGION" column
filtered_genes_data <- merge(mutations_filtered, gene_test, by = "REGION")

# Display the merged data
#print(merged_data)


# Define the desired order of categories
desired_order <- c('chrM', 'chrXVI', 'chrXV', 'chrXIV', 'chrXIII', 'chrXII', 'chrXI', 'chrX', 'chrIX', 'chrVIII', 'chrVII', 'chrVI', 'chrV', 'chrIV', 'chrIII', 'chrII', 'chrI')


# Convert category to a factor with the desired order
chromosometest$CHROM <- factor(chromosometest$CHROM, levels = desired_order)
filtered_genes_data$CHROM.x <- factor(filtered_genes_data$CHROM.x, levels = desired_order)
filtered_genes_data$CHROM.y <- factor(filtered_genes_data$CHROM.x, levels = desired_order)


# Sample data for chromosomes
chromosomes <- data.frame(
  chromosome = chromosometest$CHROM,
  length = chromosometest$length
)

# Sample data for genes
genes <- data.frame(
  chromosome = gene_test$CHROM,
  start = gene_test$START,
  end = gene_test$STOP,
  geneName = gene_test$GENE 
)


filtered_genes <- data.frame(
  chromosome = filtered_genes_data$CHROM.x,
  start = filtered_genes_data$START,
  end = filtered_genes_data$STOP,
  geneName = filtered_genes_data$GENE.y
)

filtered_genes_data$CHROM.x <- factor(filtered_genes_data$CHROM.x, levels = desired_order)


# Plotting
p <- ggplot() +
  geom_bar(data = chromosomes, aes(x = length, y = chromosome), stat = 'identity', fill = 'lightblue', width = 0.5) + # swapped x and y
  geom_rect(data = genes, aes(ymin = as.numeric(factor(chromosome)) - 0.4, # swapped ymin and ymax
                              ymax = as.numeric(factor(chromosome)) + 0.4,
                              xmin = start, 
                              xmax = end,
                              text = geneName), 
            fill = 'red', alpha = 0.5) +
  labs(title = 'Chromosomes with Genes Overlay',
       y = 'Chromosome', # changed x-axis label to Chromosome
       x = 'Length') # changed y-axis label to Length

# Convert ggplot2 plot to plotly
ggplotly(p)
