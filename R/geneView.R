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
    
    # Keeps last choice so it doesn't change as much when switching selections
    last_choices <- reactiveVal(character(0))
    
    # Updates gene dropdown only when choices actually change
    stable_filtered <- debounce(filtered_data, 200) # adding a lag so that filtered_data can actually finish updating and it doesn't grab an intermediate step
    
    observeEvent(stable_filtered(), {
      new_choices <- tryCatch({
        filtered_data() %>%
          dplyr::pull(GENE) %>%
          intersect(genes_info$GENE) %>%
          unique() %>%
          as.character() %>%
          sort()
      }, error = function(e) {
        character(0)
      })
      
      # Compare to previous choices, only updates when different
      if (!identical(new_choices, last_choices())) {
        # preserve current selection if possible
        current_sel <- isolate(input$geneSelectDropDown)
        selected_gene <- if (!is.null(current_sel) && nzchar(current_sel) && current_sel %in% new_choices) {
          current_sel
        } else if (length(new_choices) > 0) {
          new_choices[1]
        } else {
          NULL
        }
        
        updateSelectizeInput(
          session,
          "geneSelectDropDown",
          choices = as.character(new_choices),
          selected = selected_gene,
          server = TRUE,
          options = list(maxOptions = length(new_choices))
        )
        
        # store for next comparison
        last_choices(new_choices)
      }
    }, ignoreNULL = FALSE, ignoreInit = FALSE)
    
    # reactive wrapper for selected gene (keeps normal immediate reactivity)
    gene <- reactive(input$geneSelectDropDown)
    
    # Learn about Gene button within gene viewer
    if (link != "NONE") {
      output$url <- renderUI({
        req(gene())
        a("Learn about Gene",
          href = gene_info_link_function(genes_info, gene()),
          class = "btn btn-default", target = "_blank"
        )
      })
    }
    
    # Plot
    output$geneViewPlot <- renderPlotly({
      # require a gene selection and filtered data
      req(gene(), stable_filtered())
      
      ranges <- reactiveValues(x = NULL, y = NULL)
      
      # message when nothing selected
      geneview_message <- "Please select a gene below"
      gene_message_length <- nchar(geneview_message)
      gene_spaces_on_each_side <- floor((total_spaces - gene_message_length) / 2)
      select_gene_message <- sprintf(
        "%s%s%s%s",
        "\n\n\n\n\n\n\n\n",
        paste(rep(" ", gene_spaces_on_each_side), collapse = ""),
        geneview_message,
        paste(rep(" ", gene_spaces_on_each_side), collapse = "")
      )
      
      validate(
        need(gene(), select_gene_message),
        need(stable_filtered(), "Loading data...")
      )
      
      # Use a local snapshot of stable_filtered() to avoid extra invalidations mid-render
      mutation_data_value <- isolate(stable_filtered())
      
      # Merge only if both exist
      common_cols <- intersect(colnames(mutation_data_value), colnames(genes_info))
      if (length(common_cols) == 0) {
        validate("No matching columns between mutation data and genes_info.")
      }
      mutation_data_value <- merge(mutation_data_value, genes_info, by = common_cols)
      mutation_data_value <- mutation_data_value[order(mutation_data_value$GENE), ]
      
      # Filter for the gene
      cur_gene <- dplyr::filter(mutation_data_value, GENE == gene())
      
      # pattern and annotations (unchanged)
      pattern <- "(?<=\\d)([A-Za-z]|\\*|indel)$|([A-Za-z]|\\*)$"
      all_annotations <- c("missense", "nonsense", "5'-upstream", "indel-frameshift", "indel-inframe", "synonymous")
      
      # compute counts
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
          .groups = "drop"
        )
      
      if (nrow(count_proteins) == 0) {
        validate("Loading data...")
      }
      
      count_proteins_same <- count_proteins %>%
        group_by(Numbers) %>%
        summarize(
          GENE = first(GENE),
          PROTEIN = list(PROTEIN),
          ANNOTATION = first(ANNOTATION),
          Counts_diff_mutation = list(COUNTS),
          Counts_tot = sum(COUNTS),
          indel = first(indel),
          .groups = "drop"
        )
      
      count_proteins_same <- count_proteins_same %>%
        mutate(
          combined = purrr::map2_chr(
            PROTEIN, Counts_diff_mutation,
            function(prot, counts) {
              prot_list <- str_split(prot, ", ") %>% unlist()
              counts_list <- str_split(counts, ", ") %>% unlist()
              combined_strings <- purrr::map2_chr(prot_list, counts_list, function(p, c) {
                letter1 <- substr(p, 1, 1)
                letter2 <- str_extract(p, pattern)
                paste("Count", letter1, "->", letter2, ":", c, "\n", sep = " ")
              })
              paste(combined_strings, collapse = "")
            }
          )
        ) %>%
        mutate(
          PROTEIN = as.character(PROTEIN),
          AA_WT = substr(PROTEIN, 1, 1),
          AA_POS = if_else(ANNOTATION == "5'-upstream", -15, as.numeric(str_extract(PROTEIN, "[0-9]+"))),
          AA_M = substr(PROTEIN, nchar(PROTEIN), nchar(PROTEIN)),
          ANNOTATION = factor(ANNOTATION, levels = all_annotations)
        )
      
      annotation_colors <- set_names(color_vector, all_annotations)
      
      # ranges
      xmax <- genes_info %>%
        filter(GENE == gene()) %>%
        pull(PROTEIN_LENGTH) %>%
        unique() %>%
        as.numeric()
      if (length(xmax) == 0 || is.na(xmax)) xmax <- 0
      ranges$x <- c(-50, xmax + 50)
      
      ymax_count <- if (nrow(count_proteins_same) == 0) 0 else max(count_proteins_same$Counts_tot, na.rm = TRUE)
      ranges$y <- c(0, ymax_count + ifelse(ymax_count < 5, 4, 1))
      
      # build ggplot (kept your original appearance + minor safe guards)
      p <- count_proteins_same %>%
        ggplot(aes(x = AA_POS, y = Counts_tot)) +
        geom_point(data = data.frame(
          PROTEIN = rep(NA_character_, length(all_annotations)),
          AA_WT = rep(NA_character_, length(all_annotations)),
          AA_POS = rep(NA_real_, length(all_annotations)),
          AA_M = rep(NA_character_, length(all_annotations)),
          ANNOTATION = factor(all_annotations, levels = all_annotations),
          Counts_diff_mutation = I(rep(list(numeric(0)), length(all_annotations))),
          Counts_tot = rep(NA_integer_, length(all_annotations)),
          combined = rep(NA_character_, length(all_annotations))
        ), aes(y = 0, color = ANNOTATION, text = NULL), size = 2) +
        geom_hline(yintercept = 0, linetype = 2, alpha = .2) +
        geom_segment(aes(x = 0, xend = xmax, y = 0, yend = 0), size = 15, color = "cornflowerblue") +
        geom_segment(aes(x = AA_POS, xend = AA_POS, y = 0, yend = Counts_tot), color = "gray") +
        geom_point(aes(
          x = AA_POS,
          y = Counts_tot,
          color = ANNOTATION,
          text = ifelse(
            grepl("indel", ANNOTATION),
            paste0("Indel\n", abs(indel), " base ", ifelse(indel > 0, "insertion", "deletion"), "\nCount: ", Counts_diff_mutation, "\nPosition: ", AA_POS),
            ifelse(is.na(PROTEIN),
                   paste0(ANNOTATION, "\nCount: ", Counts_diff_mutation, "\nPosition: ", -abs(unique(cur_gene$POS) - unique(cur_gene$START))),
                   paste0(combined, "\nPosition: ", AA_POS)
            )
          )
        ), size = 2) +
        ggtitle(as.character(gene())) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(20, 20, 0, 20)) +
        xlab("Amino acid position") +
        ylab("Mutation Count") +
        scale_y_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
        scale_color_manual(values = annotation_colors) +
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
        annotate("text", x = 1, y = Inf, label = "Drag over mutations to see more", hjust = 0, vjust = 2, color = "black", size = 5) +
        guides(color = guide_legend(title = "Annotation"))
      
      # convert to plotly, register events so other listeners work
      ggplotly(p, tooltip = "text") %>%
        event_register('plotly_hover') %>%
        event_register('plotly_click') %>%
        event_register('plotly_relayout') %>%
        config(toImageButtonOptions = list(format = 'svg', filename = 'gene_view_plot', height = 500, width = 700))
    })
  })
}
