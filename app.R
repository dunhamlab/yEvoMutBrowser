library(shiny)
library(DBI)
library(RSQLite)

# Connect to SQLite database
con <- dbConnect(RSQLite::SQLite(), "GenomicsDatabase.db")

dbListTables(con)

# Disconnect from database

dbDisconnect(con)
#Boilerplate code
ui <- pageWithSidebar(
  
  # App title ----
  headerPanel("Please insert your files into the form below"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    fileInput("file1", "Upload VCF", accept = ".gz"),
    fileInput("file2", "Upload GFF", accept = ".gff"),
    fileInput("file3", "Upload FAStA", accept = ".fasta"),
    actionButton("action", label = "Write to DB"),
    hr()
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