PATH_TO_VCF_CSV <- "all_yEvo_vcf.csv"

ORGANISM_GENE_INFO_PATH <- "gene_info.csv"
ORGANISM_CHROMOSOME_INFO_PATH <- "chromosome_info.csv"

ORGANISM_GENE_INFO_LINK <- "https://www.yeastgenome.org/locus/"
ORGANISM_GENE_INFO_LINK_FUNCTION <- function (genes_info, selected_gene) {
  sgdid <- reactiveValues(value = NULL)
  sgdid_values <- genes_info[genes_info$GENE == selected_gene, "SGDID"]
  sgdid$value <- sgdid_values
  paste0(ORGANISM_GENE_INFO_LINK, sgdid$value)
}
