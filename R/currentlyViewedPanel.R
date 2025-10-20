current_view_panel_ui <- function(id) {
  textOutput(NS(id, "currentlyViewing"))
}


current_view_panel_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$currentlyViewing <- renderText({"Currently viewing..."})
  })
}
