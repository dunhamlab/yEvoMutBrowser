library(shiny)
library(DBI)
library(RSQLite)
library(ggplot2)
library(dplyr)

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
    plotOutput("plot")
    
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
  table <- read.table("MasterVCF.txt", header=TRUE)
  output$plot <- renderPlot({
    blank_theme <- theme_minimal()+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold")
      )
    table %>% count(CHROM) %>% mutate(percent=n/sum(n)*100) %>% ggplot(aes(x="",y=percent,fill=CHROM)) + geom_bar(stat="identity", width=1) +
      coord_polar("y", start=0) + blank_theme + ggtitle("Percentage of Variants by Chromosome")
  })
  
}

shinyApp(ui, server)