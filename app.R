library(shiny)
library(DBI)
library(RSQLite)
library(ggplot2)
library(dplyr)
library(PLColors)

## Only run examples in interactive R sessions
if (interactive()) {
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        fileInput("file1", "Choose VCF File", accept = ".csv"),
        checkboxInput("header", "Header", TRUE)
      ),
      mainPanel(
        tableOutput("contents"),
        plotOutput("plot", click = "plot_click"),
        verbatimTextOutput("info")
      )
    )
  )
  
  server <- function(input, output) {
    output$contents <- renderTable({
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "txt", "Please upload a VCF file"))
      
      table <- read.table(file$datapath, header = input$header)
      head(table)
    })
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
      table %>% count(ANNOTATION) %>% mutate(percent=n/sum(n)*100) %>% ggplot(aes(x="",y=percent,fill=ANNOTATION)) + geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) + blank_theme + ggtitle("Percentage of Variants by Type") + geom_text(aes(label = round(percent), digits = 0),position = position_stack(vjust = 0.5),col="white") +blank_theme + scale_fill_manual(values=pl_palette("lorax",5))
    })
    
    output$info <- renderPrint({
      req(input$plot_click)
      x <- round(input$plot_click$x, 2)
      y <- round(input$plot_click$y, 2)
      cat("[", x, ", ", y, "]", sep = "")
    })
  }
  
  shinyApp(ui, server)
}

shinyApp(ui, server)