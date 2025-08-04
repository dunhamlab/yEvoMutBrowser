# yEvo Mutation Browser Configuration

# DEFAULT DATA ---------------------------------------------------------------------------------------------------------

PATH_TO_VCF_CSV <- "all_yEvo_vcf.csv"

# ORGANISM--------------------------------------------------------------------------------------------------------------

ORGANISM_GENE_INFO_PATH <- "new_gene_info.csv"
ORGANISM_CHROMOSOME_INFO_PATH <- "chromosome_info.csv"
ORGANISM_GENE_INFO_LINK <- "https://www.yeastgenome.org/locus/"
ORGANISM_GENE_INFO_LINK_FUNCTION <- function (genes_info, selected_gene) {
  sgdid <- reactiveValues(value = NULL)
  sgdid_values <- genes_info[genes_info$GENE == selected_gene, "SGDID"]
  sgdid$value <- sgdid_values
  paste0(ORGANISM_GENE_INFO_LINK, sgdid$value)
}

# COLORS----------------------------------------------------------------------------------------------------------------

VARIANTS_PIE_CHART_COLORS <- c(
  "#9edae5", "#17becf", "#dbdb8d", "#bcbd22", "#c7c7c7",
  "#e377c2", "#7f7f7f", "#f7b6d2", "#c49c94", "#8c564b",
  "#c5b0d5", "#9467bd", "#ff9896", "#d62728", "#98df8a",
  "#2ca02c", "#ffbb78", "#ff7f0e", "#aec7e8", "#1f77b4"
)
SNP_CHART_COLORS <- c("#0072B2", "#CC79A7", "#009E73")
GENE_VIEW_COLORS <- c("#800080", "#61D04F", "#F5C710", "#DF536B", "#FFB6C1", "cornflowerblue")