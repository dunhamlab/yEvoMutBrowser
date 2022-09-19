library(shiny)

#Boilerplate code
ui <- pageWithSidebar(
  
  # App title ----
  headerPanel("Miles Per Gallon"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    
    #Input goes here
  ),
  
  
  # Main panel for displaying outputs ----
  mainPanel(
    
    #Output goes here
  )
)

#Preprocess data here

server <- function(input, output) {
  #Put the actual logic for vcf logic and visualizations here
}

shinyApp(ui, server)