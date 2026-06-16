protein_prediction_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
           "Protein Prediction",
           div(style = "margin-left: 20px; width: 600px;",
           div("", style = "height: 10px;"),
           tags$div(class = "info_div",
             verbatimTextOutput(NS(id, "info")),
           ),
            
            selectInput(NS(id, "geneSelectDropDown"), "Gene", choices = NULL),

            # IGV-style sequence track: a canvas drawn by static/sequence-track.js.
            # Scroll to zoom, drag to pan; genomic coordinates show on hover.
            tags$div(
              id = NS(id, "seqtrack-wrap"),
              style = paste(
                "position: relative; width: 100%; height: 95px;",
                "border: 1px solid #ddd; border-radius: 4px; overflow: hidden;",
                "user-select: none;"
              ),
              tags$canvas(
                id = NS(id, "seqtrack-canvas"),
                style = "width: 100%; height: 100%; display: block;"
              ),
              tags$div(
                id = NS(id, "seqtrack-tip"),
                style = paste(
                  "position: absolute; display: none; pointer-events: none;",
                  "background: rgba(0,0,0,0.8); color: #fff; padding: 2px 6px;",
                  "border-radius: 3px; font: 12px monospace; white-space: nowrap;",
                  "z-index: 10;"
                )
              )
            ),

            tags$script(src = "static/sequence-track.js"),

           ))
}

protein_prediction_server <- function(id, total_spaces, filtered_data, genes_info, link, gene_info_link_function, color_vector, alphafold_color_vector) {
  moduleServer(id, function(input, output, session) {
    output$info <- renderText({
      "This section provides a visual representation of the predicted protein structure based on the gene sequence. The reference sequence is derived from S288c."
    })

    stable_data <- debounce(filtered_data, 150)
    last_choices <- reactiveVal(NULL)

    #generates information based on the stable data (so that it doesn't flip back and forth during filtering)
    observeEvent(stable_data(), {

      new_choices <- tryCatch({
        stable_data() %>%
          dplyr::pull(GENE) %>%
          intersect(genes_info$GENE) %>%
          unique() %>%
          sort() %>%
          as.character()
      }, error = function(e) character(0))

      if (identical(new_choices, last_choices()))
        return()

      #holds onto selected gene
      current_sel <- isolate(input$geneSelectDropDown)

      selected_gene <- if (!is.null(current_sel) &&
                             nzchar(current_sel) &&
                             current_sel %in% new_choices) {
        current_sel
      } else if (length(new_choices) > 0) {
        new_choices[1]
      } else {
        NULL
      }

      updateSelectizeInput(
        session,
        "geneSelectDropDown",
        choices = new_choices,
        selected = selected_gene,
        server = TRUE,
        options = list(maxOptions = length(new_choices))
      )

      last_choices(new_choices)

    }, ignoreInit = FALSE)

    gene_sequence <- read.csv("data/gene_sequences.csv")

    gene <- reactive({input$geneSelectDropDown})

    cur_gene <- reactive({
      # merge & filter gene info and mutation info and filter for selected gene, using the stablized data
      fd <- stable_data() %>%
        merge(genes_info, by = intersect(names(stable_data()), names(genes_info))) %>%
        arrange(GENE)

      filter(fd, GENE == gene())
    })
    # observeEvent(gene(), {
    #   print(cur_gene())
    # })

    # Push the selected gene's reference sequence to the canvas renderer
    # (static/sequence-track.js). The JS handles all drawing, zoom and pan,
    # so only the raw sequence + genomic anchor cross the wire.
    observeEvent(gene(), {
      req(gene())
      ids <- list(
        canvasId = session$ns("seqtrack-canvas"),
        tipId    = session$ns("seqtrack-tip")
      )
      matched_row <- gene_sequence[gene_sequence$GENE == gene(), ]

      if (nrow(matched_row) < 1) {
        session$sendCustomMessage("seqTrackData", c(ids, list(
          empty = TRUE,
          message = paste0("No reference sequence available for ", gene(), ".")
        )))
        return()
      }

      matched_row <- matched_row[1, ]  # guard against duplicate rows
      seq_str <- toupper(as.character(matched_row$sequence))

      if (is.na(seq_str) || !nzchar(seq_str)) {
        session$sendCustomMessage("seqTrackData", c(ids, list(
          empty = TRUE,
          message = paste0("Sequence for ", gene(), " is empty.")
        )))
        return()
      }

      session$sendCustomMessage("seqTrackData", c(ids, list(
        seq   = seq_str,
        start = as.numeric(matched_row$start),
        chrom = as.character(matched_row$chrom)
      )))
    }, ignoreInit = FALSE)
})
}