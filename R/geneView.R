gene_view_ui <- function(id) {
  tabPanel(
    "Gene View", div("", style = "height: 10px;"),
    plotlyOutput(NS(id, "geneViewPlot"), width = "600px"), verbatimTextOutput(NS(id, "gene")),
    selectInput(NS(id, "geneSelectDropDown"), "Gene", choices = NULL),

    uiOutput(NS(id, "url"))
  )
}

gene_view_server <- function(id, total_spaces, filtered_data, genes_info, link, gene_info_link_function, color_vector) {
  moduleServer(id, function(input, output, session) {
    #gene view draopdown menu
    observe({
      choices <- filtered_data() %>% pull(GENE)

      # Filter choices to include only those present in genes_info
      choices <- choices[choices %in% genes_info$GENE]

      choices <- sort(choices)

      updateSelectizeInput(session, "geneSelectDropDown",
                           choices = as.character(choices),
                           server = TRUE,
                           options = list(maxOptions = length(choices))
      )
    })

    gene <- reactive(input$geneSelectDropDown)

    # Learn about Gene button within gene viewer
    if(link != "NONE") {
    output$url <- renderUI({
      url <- a("Learn about Gene",
               href = gene_info_link_function(genes_info, input$geneSelectDropDown),
               class = "btn btn-default", target = "_blank"
      )
      url
    })
    }


    output$geneViewPlot <- renderPlotly({
      ranges <- reactiveValues(x = NULL, y = NULL)

      # Showing please select gene message
      # Construct the string with spaces on each side of the loading message
      geneview_message <- "Please select a gene below"
      # Calculate the number of empty spaces needed on each side
      gene_message_length <- nchar(geneview_message)
      gene_spaces_on_each_side <- floor(
        (total_spaces - gene_message_length) / 2
      )
      select_gene_message <- sprintf(
        "%s%s%s%s",
        "\n\n\n\n\n\n\n\n",
        paste(rep(" ", gene_spaces_on_each_side), collapse = ""),
        geneview_message,
        paste(rep(" ", gene_spaces_on_each_side), collapse = "")
      )

      # Render gene view plot

      validate(
        need(gene(), select_gene_message),
        need(filtered_data(), "Loading data...")
      )

      mutation_data_value <- filtered_data()

      # Merge the data frames based on the common columns
      common_cols <- intersect(
        colnames(mutation_data_value),
        colnames(genes_info)
      )
      mutation_data_value <- merge(mutation_data_value, genes_info,
                                   by = common_cols
      )
      mutation_data_value <- mutation_data_value[order(
        mutation_data_value$GENE
      ),]

      # Filter data for the specific gene
      cur_gene <- filter(mutation_data_value, GENE == gene())

      # Pattern for extracting the second part of protein
      pattern <- "(?<=\\d)([A-Za-z]|\\*|indel)$|([A-Za-z]|\\*)$"

      all_annotations <- c(
        "missense", "nonsense", "5'-upstream",
        "indel-frameshift", "indel-inframe", "synonymous", "transposon"
      )
      # IF ADDING NEW ANNOTATIONS - DON'T FORGET TO ADD BOTH HERE AND IN
      # annotation_colors BELOW

      # Group and summarize protein counts
      count_proteins <- cur_gene %>%
        mutate(indel = nchar(ALT) - nchar(REF)) %>% # Calculate indel difference
        group_by(POS, PROTEIN) %>%
        summarize(
          GENE = first(GENE),
          PROTEIN = first(PROTEIN),
          ANNOTATION = first(ANNOTATION),
          COUNTS = n(),
          Letter1 = substr(PROTEIN, 1, 1), # Extract the first character
          # Extract the numbers
          Numbers = as.numeric(str_extract(PROTEIN, "[0-9]+")),
          Letter2 = str_extract(PROTEIN, pattern),
          indel = first(indel),
          POS = first(POS),
          START = first(START)
        ) %>%
        ungroup()

      # If count_proteins$COUNT is empty, the data is not fully loaded in yet
      if (length(count_proteins$COUNTS) <= 0) {
        validate("Loading data...")
      }

      # Group by position numbers and summarize
      count_proteins_same <- count_proteins %>%
        group_by(Numbers) %>%
        summarize(
          GENE = first(GENE),
          PROTEIN = list(PROTEIN),
          ANNOTATION = first(ANNOTATION),
          Counts_diff_mutation = list(COUNTS),
          Counts_tot = sum(COUNTS),
          indel = first(indel),
          POS = first(POS),
          START = first(START)
        ) %>%
        ungroup()

      # Combine protein and count strings
      count_proteins_same <- count_proteins_same %>%
        mutate(
          combined = map2_chr(
            PROTEIN, Counts_diff_mutation,
            function(prot, counts) {
              prot_list <- str_split(prot, ", ") %>% unlist()
              counts_list <- str_split(counts, ", ") %>% unlist()

              combined_strings <- map2_chr(
                prot_list, counts_list,
                function(p, c) {
                  letter1 <- substr(p, 1, 1)
                  # numbers <- str_extract(p, "[0-9]+") %>% as.numeric()
                  letter2 <- str_extract(p, pattern)
                  paste("Count", letter1, "->", letter2, ":", c, "\n",
                        sep = " "
                  )
                }
              )
              paste(combined_strings, collapse = "")
            }
          )
        ) %>%
        mutate(
          PROTEIN = as.character(PROTEIN),
          # Extract the first character Amino Acid Wild Type
          AA_WT = substr(PROTEIN, 1, 1),
          AA_POS = if_else(ANNOTATION == "5'-upstream", -15,
                           if_else(ANNOTATION == "transposon",
                                   # Calculate amino acid position from nucleotide position for transposons
                                   as.numeric(ceiling((POS - START + 1) / 3)),
                                   # Extract Amino Acid Position for regular mutations
                                   as.numeric(str_extract(PROTEIN, "[0-9]+"))
                           )
          ),
          # Amino Acid Mutation
          AA_M = substr(PROTEIN, nchar(PROTEIN), nchar(PROTEIN)),
          ANNOTATION = factor(ANNOTATION, levels = all_annotations)
        )

      annotation_colors <- set_names(color_vector, all_annotations)

      # Generating ranges for the plot size
      xmax <- genes_info %>%
        filter(GENE == gene()) %>%
        pull(PROTEIN_LENGTH) %>%
        unique() %>%
        as.numeric()
      ranges$x <- c(-50, xmax + 50)

      ranges$y <- c(0, max(count_proteins_same$Counts_tot) +
        ifelse(max(count_proteins_same$Counts_tot) < 5, 4, 1))

      p <- count_proteins_same %>%
        ggplot(
          aes(
            x = AA_POS,
            y = Counts_tot,
          )
        ) +
        # dummy geom_point to get legend showing all_annotations regardless of
        # what data is displayed
        geom_point(data = data.frame(
          PROTEIN = rep(NA_character_, length(all_annotations)),
          AA_WT = rep(NA_character_, length(all_annotations)),
          AA_POS = rep(NA_real_, length(all_annotations)),
          AA_M = rep(NA_character_, length(all_annotations)),
          ANNOTATION = factor(all_annotations, levels = all_annotations),
          Counts_diff_mutation = I(rep(
            list(numeric(0)),
            length(all_annotations)
          )),
          Counts_tot = rep(NA_integer_, length(all_annotations)),
          combined = rep(NA_character_, length(all_annotations))
        ), aes(y = 0, color = ANNOTATION, text = NULL), size = 2) +
        geom_hline(yintercept = 0, linetype = 2, alpha = .2, aes(text = NULL)) +
        geom_segment(aes(x = 0, xend = xmax, y = 0, yend = 0, text = NULL),
                     size = 15, color = "cornflowerblue"
        ) +
        geom_segment(aes(
          x = AA_POS, xend = AA_POS, y = 0, yend = Counts_tot,
          text = NULL
        ), color = "gray") +
        geom_point(aes(
          x = AA_POS,
          y = Counts_tot,
          color = ANNOTATION,
          text = ifelse(
            ANNOTATION == "transposon", # Check if ANNOTATION is transposon
            paste0(
              "Transposon insertion\nPosition: ", AA_POS
            ), # Text for transposon annotations
            ifelse(
              grepl("indel", ANNOTATION), # Check if ANNOTATION contains "indel"
              paste0(
                "Indel\n", abs(indel), " base ", ifelse(indel > 0,
                                                        "insertion", "deletion"
                ), "\nCount: ", Counts_diff_mutation,
                "\nPosition: ", AA_POS
              ), # Text for indel annotations
              ifelse(
                is.na(PROTEIN),
                paste0(
                  ANNOTATION, "\nCount: ", Counts_diff_mutation,
                  "\nPosition: ", -abs(unique(cur_gene$POS) -
                                         unique(cur_gene$START))
                ),
                paste0(combined, "\nPosition: ", AA_POS)
              )
            )
          )
        ), size = 2) +
        ggtitle(as.character(gene())) +
        theme_classic() +
        theme(
          axis.text.x = element_text(),
          axis.title.y = element_text(),
          axis.text.y = element_text(),
          axis.ticks.y = element_line(),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(20, 20, 0, 20)
        ) +
        xlab("Amino acid position") +
        ylab("Mutation Count") +
        scale_y_continuous(breaks = function(x) {
          seq(floor(min(x)),
              ceiling(max(x)),
              by = 1
          )
        }) +
        scale_color_manual(values = annotation_colors) + # Manual color scale
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
        annotate("text",
                 x = 1, y = Inf,
                 label = "Drag over mutations to see more", hjust = 0, vjust = 2,
                 color = "black", size = 5
        ) +
        guides(color = guide_legend(title = "Annotation"))

      ggplotly(p, tooltip = "text") %>%
        config(
          toImageButtonOptions = list(
            format = 'svg',
            filename = 'gene_view_plot',
            height = 500,
            width = 700
          )
        )
    })
  })
}