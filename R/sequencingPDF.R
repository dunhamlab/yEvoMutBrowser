sequencing_pdf_ui <- function(id) {
  tabPanel(
    "How Does Sequencing Work?",
    uiOutput(NS(id, "pdf_viewer"))
  )
}

sequencing_pdf_server <- function(id, pdf_path) {
  moduleServer(id, function(input, output, session) {
    output$pdf_viewer <- renderUI({
      tags$iframe(
        style = "height:1000px;width:100%;scrolling:yes",
        src = pdf_path
      )
    })
  })
}