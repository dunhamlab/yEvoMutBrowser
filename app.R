library(shiny)
library(shinythemes)
library(DBI)
library(RSQLite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(PLColors)
library(forcats)

table <- read.table("final_MASTERVCF.txt", header=TRUE)
addResourcePath(prefix = 'img', directoryPath = '~/GSHackathon/img')
#addResourcePath("~/GSHackathon", "GSHackathon")
## Only run examples in interactive R sessions
if (interactive()) {

  ui <-  navbarPage(
                    title = div(img(src="img/yEvo_logo.png",
                                    filetype = "image/png",
                                    style="margin-top: -14px;
                               padding-right:10px;
                               padding-bottom:10px",
                                    height = 60)
                                ),
                    #theme = "journal",
                    windowTitle="yEvo",
                    theme = shinytheme("cerulean"),
                    
                    
                    
                   tabPanel("Data Visualizations",
    #fluidPage(
    sidebarLayout(
      sidebarPanel(
        fileInput("file1", "Choose VCF File", accept = ".txt"),
        checkboxInput("header", "Header", TRUE),
        selectInput("SAMPLE", "Sample",
                    choices = c("RBcaff_1blue1_S9", 
                                "RBcaff1gray1_S12", 
                                "RBcaff_1orange1_S10",
                                "RBcaff_1pink2_S7",
                                "RBcaff_1yellow1_S8",
                                "RBcaff_3blue1_S15",
                                "RBcaff_3gray1_S18",
                                "RBcaff_3orange3_S16",
                                "RBcaff_3pink3_S13",
                                "RBcaff_3purple1_S17",
                                "RBcaff_3yellow1_S14",
                                "RBcaff_5blue1_S21",
                                "RBcaff_5gray3_S24",
                                "RBcaff_5orange2_S22",
                                "RBcaff_5pink2_S19",
                                "RBcaff_5purple3_S23",
                                "RBcaff_5yellow3_S20",
                                "RBcaff_6blue1_S27",
                                "RBcaff_6gray3_S30",
                                "RBcaff_6orange1_S28",
                                "RBcaff_6pink1_S25",
                                "RBcaff_6purple1_S29",
                                "RBcaff_6yellow1_S26")),
        br(),
       # actionButton("ChromosomeMap", "Chromosome Map"),
       # br(),
      #  actionButton("PieChart", "Pie Chart"),
      #  br(),
      #  actionButton("SNP", "SNP counts")
        
      ),
      
      mainPanel(
        #tableOutput("contents"),
        #verbatimTextOutput("info"),
        # Output: Tabset w/ plot, summary, and table ----
        tabsetPanel(type = "tabs",
                    tabPanel("Table", tableOutput("contents")),
                    tabPanel("Chromosome Map", plotOutput("plot1", click = "plot_click")),
                    tabPanel("Variant Pie Chart", plotOutput("plot", click = "plot_click")),
                    tabPanel("SNP Counts", plotOutput("plot2", click = "plot_click"))
        ),
        #plotOutput("plot1", click = "plot_click"),
        #plotOutput("plot", click = "plot_click"),
        #plotOutput("plot2", click = "plot_click"),

      )
    )
  ),
  tabPanel("Background",
           uiOutput("pdf_viewer") ),
  tabPanel("Page3")
                   )
  server <- function(input, output) {
    
    output$pdf_viewer <- renderUI({ tags$iframe(
      style="height:1000px;width:100%;scrolling=yes",
      src = "Black_box.pdf") }) 
    
    output$contents <- renderTable({
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "txt", "Please upload a VCF file"))
      
      table <- read.table(file$datapath, header = input$header)

      table %>% filter(SAMPLE==input$SAMPLE[1])
    })
    
    #observeEvent(input$ChromosomeMap,{
    #  cat("Showing")
    #})
    
   # observeEvent(input$PieChart,{
   #   cat("Showing")
  #  })
    
   # observeEvent(input$SNP,{
  #    cat("Showing")
   # })
    
    
    df <- eventReactive(input$ChromosomeMap,{
      table %>% 
        mutate(Chromosome=forcats::fct_relevel(Chromosome,'I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','M')) %>%
        ggplot() +
        facet_grid(vars(Chromosome),switch = "y") +
        geom_segment(aes(color=Chromosome),x = 1, y = 0, xend = table$Length, yend = 0, size=4.1,lineend = "round") +
        geom_segment(x = table$cent1, y = 0, xend = table$cent2, yend = 0, size=4.1,lineend = "round", color="black") +
        scale_color_manual(values=pl_palette("lorax",17)) +
        geom_point(aes(x=POS,y=0),shape = "|", size=2.9, data=table
                   %>% mutate(Chromosome=forcats::fct_relevel(Chromosome,'I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','M'))%>% 
                     filter(SAMPLE==input$SAMPLE[1])) + 
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.background = element_blank()
        ) +
        xlim(c(0,1540000)) +
        ggtitle("Where do variants fall on chromosomes") + 
        xlab("Position along chromosome") + ylab("Chromosome") +
        theme(strip.text.y.left = element_text(angle = 0),
              plot.title = element_text(hjust = 0.5),
              legend.position="none")
    })
    
    df1 <- eventReactive(input$PieChart,{
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

      table %>% filter(SAMPLE==input$SAMPLE[1]) %>% 
        count(ANNOTATION) %>% 
        mutate(percent=n/sum(n)*100) %>% 
        ggplot(aes(x="",y=percent,fill=ANNOTATION)) + 
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) + 
        ggtitle("Percentage of Variants by Type") + 
        geom_text(aes(label = round(percent), digits = 0),position = position_stack(vjust = 0.5),col="white") + 
        blank_theme + 
        scale_fill_manual(values=pl_palette("lorax",5))
    })
    
    df2 <- eventReactive(input$SNP,{
      
      table %>% mutate(transition=paste(REF,"_",ALT, sep=""))  %>% 
        select(transition,SAMPLE) %>% 
        filter(SAMPLE==input$SAMPLE[1]) %>%
        mutate(length = nchar(transition)) %>% 
        filter(length == 3) %>%  
        ggplot(aes(x=as.factor(transition),fill=as.factor(transition))) + 
        geom_bar() + theme_bw() + 
        scale_fill_manual(values=pl_palette("lorax",12)) + 
        theme(legend.position = "none",strip.text.y.left = element_text(angle = 0),plot.title = element_text(hjust = 0.5)) + 
        ggtitle("Single Nucleotide transitions")  + xlab("SNP call")
      })
    

    
    output$plot1 <- renderPlot({
      #df()
      table %>% 
        mutate(Chromosome=forcats::fct_relevel(Chromosome,'I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','M')) %>%
        ggplot() +
        facet_grid(vars(Chromosome),switch = "y") +
        geom_segment(aes(color=Chromosome),x = 1, y = 0, xend = table$Length, yend = 0, size=4.1,lineend = "round") +
        geom_segment(x = table$cent1, y = 0, xend = table$cent2, yend = 0, size=4.1,lineend = "round", color="black") +
        scale_color_manual(values=pl_palette("lorax",17)) +
        geom_point(aes(x=POS,y=0),shape = "|", size=2.9, data=table
                   %>% mutate(Chromosome=forcats::fct_relevel(Chromosome,'I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','M'))%>% 
                     filter(SAMPLE==input$SAMPLE[1])) + 
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.background = element_blank()
        ) +
        xlim(c(0,1540000)) +
        ggtitle("Where do variants fall on chromosomes") + 
        xlab("Position along chromosome") + ylab("Chromosome") +
        theme(strip.text.y.left = element_text(angle = 0),
              plot.title = element_text(hjust = 0.5),
              legend.position="none")
      })
    output$plot <- renderPlot({
      #df1()
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
      
      table %>% filter(SAMPLE==input$SAMPLE[1]) %>% 
        count(ANNOTATION) %>% 
        mutate(percent=n/sum(n)*100) %>% 
        ggplot(aes(x="",y=percent,fill=ANNOTATION)) + 
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) + 
        ggtitle("Percentage of Variants by Type") + 
        geom_text(aes(label = round(percent), digits = 0),position = position_stack(vjust = 0.5),col="white") + 
        blank_theme + 
        scale_fill_manual(values=pl_palette("lorax",5))
    })
    output$plot2 <- renderPlot({
      #df2()
      
      num <- final %>% mutate(transition=paste(REF,"_",ALT, sep=""))  %>% select(transition,SAMPLE) %>% mutate(length = nchar(transition)) %>% 
        filter(SAMPLE==input$SAMPLE[1]) %>%
        count(transition) %>% 
        summarise(n = n()) %>% as.numeric()
      
      
      table %>% mutate(transition=paste(REF,"_",ALT, sep=""))  %>% 
        select(transition,SAMPLE) %>% 
        filter(SAMPLE==input$SAMPLE[1]) %>%
        mutate(length = nchar(transition)) %>% 
        filter(length == 3) %>%  
        ggplot(aes(x=as.factor(transition),fill=as.factor(transition))) + 
        geom_bar() + theme_bw() + 
        scale_fill_manual(values=pl_palette("lorax",num)) + 
        theme(legend.position = "none",
              strip.text.y.left = element_text(angle = 0),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) + 
        ggtitle("Single Nucleotide transitions")  + xlab("SNP call")
    })
  }
  
  shinyApp(ui, server)
}
