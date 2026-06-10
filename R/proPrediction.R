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

    gene <- reactive({input$geneSelectDropDown})

})}