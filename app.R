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
    dateInput("date1", "Current Date"),
    textInput("value1", "Party of Interest"),
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
#Simplified connection
server <- function(input, output) {
  #Put the actual logic for vcf logic and visualizations here
  data <- eventReactive(input$action, {
    con <- dbConnect(SQLite(), dbname="GenomicsDatabase.db")
    dbWriteTable(con, "VCF",savemode = "u", data.frame(value1 = input$file1, 
                                                          value2 = input$date1, 
                                         value3 = input$value1, 
                                         stringsAsFactors = FALSE), append = TRUE)
    data <- dbReadTable(con, "VCF")
    dbDisconnect(con)
    return(data)
  })
}

shinyApp(ui, server)