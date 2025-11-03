data_table_ui <- function(id) {
  tabPanel("Table", DT::dataTableOutput(NS(id, "data_table")))
}

data_table_server <- function(id, filtered_data) {
  moduleServer(id, function(input, output, session) {
    output$data_table <- DT::renderDataTable({
      filtered_data() %>%
        select(
          CHROM, POS, ANNOTATION, GENE, PROTEIN, condition, instructor,
          year, sample, REF, ALT
        )
    })
  })
}
