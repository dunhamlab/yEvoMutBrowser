# If you are changing the 'masterfile' that the app is built off of (aka the
# initial database),
# This needs to be changed in two locations:
# mut_backend before the server is called, and mutation_data right after the
# server is called.

### Imports
# loading necessary libraries
library(DBI)
library(devtools)
library(dplyr)
library(forcats)
library(ggplot2) ## visualization data
library(ggrepel)
library(plotly)
library(purrr)
library(RSQLite)
library(scales)
library(shiny)
library(shinyjs)
library(shinythemes)
library(stringr)
library(tidyr) ## handling data, df must be tidy to use a lot of packages
library(viridis)
library(zeallot)

PATH_TO_VCF_CSV <- "all_yEvo_vcf.csv"

yEvoMutBrowser <- function(...) {
  # loading in the VCF file to display initial choices, later turns into
  # reactive val called mutation_data that includes manually updated data
  # THIS SHOULD NOT BE CHANGED IN THE CODE. If the overall master shifts you can
  # modify it here, but mut_backend should not be written to anywere in the code
  # it is designed to just fill in in the beginning, right before the reactive
  # frame gets created.
  mut_backend <- read.csv(PATH_TO_VCF_CSV)

  # loading in the genes data file
  genes_info <- read.csv("gene_info.csv")

  # loading in the chromosomes data file
  chrom_info <- read.csv("chromosome_info.csv")

  # need to add this to upload the yEvo icon the theme
  addResourcePath(prefix = "img", directoryPath = "img")

  # create a variable called link that stores the base SGD database for locus
  link <- "https://www.yeastgenome.org/locus/"

  # how will our layout look like for the app
  ui <- navbarPage(
    useShinyjs(),
    title = div(img(
      src = "img/yEvo_logo.png",
      filetype = "image/png",
      style = "margin-top: -14px;
                          padding-right:10px;
                          padding-bottom:10px",
      height = 60
    )),
    windowTitle = "yEvo",
    theme = shinytheme("cerulean"),
    tabPanel(
      "Data Visualizations",
      sidebarLayout(
        # Left side, Class vs Cumulative View and options
        selection_panel_ui("selectionPanel", mut_backend),
        # Right side, Data Visualization
        mainPanel(
          tabsetPanel(
            type = "tabs",
            chrom_map_ui("chromMap"),
            variants_ui("variants"),
            snp_count_ui("snpCount"),
            gene_view_ui("geneView"),
            data_table_ui("dataTable"),
          )
        )
      )
    ),
    tutorial_ui("tutorial"),
    sequencing_pdf_ui("sequencingPDF")
  ) # END OF UI

  # Now entering server, which handles everything dynamically
  server <- function(input, output, session) {
    # initially setting default file of all mutation data
    mutation_data <- reactiveVal(read.csv(PATH_TO_VCF_CSV))

    shinyjs::hide("cumulDropdowns") # Initially hide cumulative drop downs

    c(selected_instructor, selected_year, selected_sample, selected_condition, selected_background, filtered_data) %<-%
      selection_panel_server("selectionPanel", filtered_data, mutation_data)

    # to create loading message below:
    loading_message <- "Loading..."
    # Calculate the number of empty spaces needed on each side
    total_spaces <- 160 # Total number of characters to occupy the line
    message_length <- nchar(loading_message)
    spaces_on_each_side <- floor((total_spaces - message_length) / 2)
    # Construct the string with spaces on each side of the loading message
    formatted_loading_message <- sprintf(
      "%s%s%s%s",
      "\n\n\n\n\n\n\n\n",
      paste(rep(" ", spaces_on_each_side), collapse = ""),
      loading_message,
      paste(rep(" ", spaces_on_each_side), collapse = "")
    )

    # Define the desired order of categories
    desired_order <- c(
      "chrM", "chrXVI", "chrXV", "chrXIV", "chrXIII",
      "chrXII", "chrXI", "chrX", "chrIX", "chrVIII", "chrVII", "chrVI", "chrV",
      "chrIV", "chrIII", "chrII", "chrI"
    )
    # Convert category to a factor with the desired order
    chrom_info$CHROM <- factor(chrom_info$CHROM, levels = desired_order)

    # Define a mapping from chromosome names to numbers
    chromosome_mapping <- c(
      chrM = 1, chrXVI = 2, chrXV = 3, chrXIV = 4, chrXIII = 5, chrXII = 6,
      chrXI = 7, chrX = 8, chrIX = 9, chrVIII = 10, chrVII = 11, chrVI = 12,
      chrV = 13, chrIV = 14, chrIII = 15, chrII = 16, chrI = 17
    )

    # Create an empty dataframe to store the final results
    final_gene <- reactive({
      mutation_data_value <- filtered_data()
      # Merge the data frames based on the “REGION” column
      common_cols <- intersect(
        colnames(mutation_data_value),
        colnames(genes_info)
      )
      mutation_data_value <- merge(mutation_data_value, genes_info,
                                   by = common_cols
      )

      # Iterate through unique genes
      final_gene_static <- mutation_data_value %>%
        group_by(GENE) %>%
        summarize(
          CHROM = first(CHROM),
          START = first(START),
          STOP = first(STOP),
          Counts = n(),
          chrom_as_num = first(chromosome_mapping[match(
            first(CHROM),
            names(chromosome_mapping)
          )])
        ) %>%
        ungroup()

      return(final_gene_static)
    })

    # Render the dataframe in the tableOutput
    data_table_server("dataTable", filtered_data)

    sequencing_pdf_server("sequencingPDF", "img/Black_box.pdf")

    chrom_map_server(
      "chromMap", final_gene,
      formatted_loading_message, chrom_info
    )
    variants_server(
      "variants", mutation_data,
      filtered_data
    )
    snp_count_server("snpCount", filtered_data)
    gene_view_server(
      "geneView", total_spaces, filtered_data, genes_info, link
    )

    observeEvent(input$append_btn, {
      new_csv_path <- input$new_csv$datapath
    })
  }

  shinyApp(ui, server)
}
