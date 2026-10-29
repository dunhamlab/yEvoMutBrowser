# yEvo Mutation Browser Configuration

# DEFAULT DATA ---------------------------------------------------------------------------------------------------------

PATH_TO_VCF_CSV <- "all_yEvo_vcf.csv"

# ORGANISM--------------------------------------------------------------------------------------------------------------

ORGANISM_GENE_INFO_PATH <- "gene_info.csv"
ORGANISM_CHROMOSOME_INFO_PATH <- "chromosome_info.csv"
ORGANISM_PFAM_DOMAIN_INFO_PATH <- "pfam.csv"
ORGANISM_PROSITE_DOMAIN_INFO_PATH <- "prosite.csv"
ORGANISM_PATHWAY_INFO_PATH <- "pathways.csv"
ORGANISM_MOTIF_INFO_PATH <- "motifs.csv"
GENE_SEQUENCE_PATH <- "data/gene_sequences.csv"

ORGANISM_GENE_INFO_LINK <- "https://www.yeastgenome.org/locus/"

ORGANISM_GENE_INFO_LINK_FUNCTION <- function (genes_info, selected_gene) {
  sgdid <- reactiveValues(value = NULL)
  sgdid_values <- genes_info[genes_info$GENE == selected_gene, "SGDID"]
  sgdid$value <- sgdid_values
  paste0(ORGANISM_GENE_INFO_LINK, sgdid$value)
}



# COLORS----------------------------------------------------------------------------------------------------------------

VARIANTS_PIE_CHART_COLORS <- c(
  "#F5C710", "#17becf", "#DF536B", "#f7b6d2",  "#c49c94",
  "#8c564b", "#800080", "#61D04F", "#7f7f7f", "cornflowerblue",
  "#c5b0d5", "#9467bd", "#ff9896", "#d62728", "#98df8a", "#e377c2",
  "#2ca02c", "#ffbb78", "#ff7f0e", "#aec7e8", "#1f77b4", "#9edae5"
)
SNP_CHART_COLORS <- c("#0072B2", "#CC79A7", "#009E73")
GENE_VIEW_COLORS <- c("#800080", "#61D04F", "#F5C710", "#DF536B", "#FFB6C1", "#6495ED", "#6A5ACD")

# Shared mutation annotation colors.
# The first seven colors preserve the existing Gene Pro View
# and Protein Prediction annotation scheme.
ANNOTATION_COLORS <- stats::setNames(
  GENE_VIEW_COLORS,
  c(
    "missense",
    "nonsense",
    "5'-upstream",
    "indel-frameshift",
    "indel-inframe",
    "synonymous",
    "transposon"
  )
)

# Transparency for mutation annotation colors:
# 0 = fully transparent, 1 = fully opaque
ANNOTATION_ALPHA <- 0.4

# Return shared colors for requested annotations,
# applying the configured transparency.
# Unknown annotation types receive deterministic additional colors.
get_annotation_colors <- function(annotation_names) {
  
  annotation_names <- unique(as.character(annotation_names))
  annotation_names <- annotation_names[
    !is.na(annotation_names) & nzchar(annotation_names)
  ]
  
  extra_annotations <- setdiff(
    annotation_names,
    names(ANNOTATION_COLORS)
  )
  
  extra_colors <- stats::setNames(
    vapply(extra_annotations, function(annotation) {
      
      chars <- utf8ToInt(annotation)
      
      # Derive a repeatable hue from the annotation name.
      hue <- sum(as.double(chars) * seq_along(chars)) %% 360
      
      grDevices::hcl(
        h = hue,
        c = 65,
        l = 60
      )
      
    }, character(1)),
    extra_annotations
  )
  
  color_map <- c(ANNOTATION_COLORS, extra_colors)
  
  # Apply transparency directly to the requested colors.
  grDevices::adjustcolor(
    color_map[annotation_names],
    alpha.f = ANNOTATION_ALPHA
  ) |>
    stats::setNames(annotation_names)
}

ALPHAFOLD_COLORS <- c("#0053D6", "#65CBF3", "#FFDB13", "#FF7D45")