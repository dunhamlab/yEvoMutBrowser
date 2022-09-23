library(devtools)
library(shiny)
library(shinythemes)
library(DBI)
library(RSQLite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(PLColors)
library(forcats)
library(ggrepel)
library(purrr)

final <- read.table("final_MASTERVCF.txt", header=TRUE)

addResourcePath(prefix = 'img', directoryPath = 'img')

link = "https://www.yeastgenome.org/locus/"


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
        
        selectInput("GENE", "Gene",
                    choices = c('')),
        
        selectInput("SGDID", "SGDID",
                    choices = c('')),
        
       
        uiOutput("url")
      ),
      
      mainPanel(

        tabsetPanel(type = "tabs",
                    tabPanel("Chromosome Map", plotOutput("plot1", brush =brushOpts(id = "plot_brush", fill = "#ccc", direction = "x")),verbatimTextOutput("info")),
                    tabPanel("Variant Pie Chart", plotOutput("plot", click = "plot_click"),verbatimTextOutput("text")),
                    tabPanel("SNP Counts", plotOutput("plot2", click = "plot_click")),
                    tabPanel("Gene View", value="Geneview",plotOutput("plot3", click = "plot_click")),
                    tabPanel("Table", tableOutput("contents"))

        ),


      )
    )
  ),
  tabPanel("Background",
           uiOutput("pdf_viewer") ),
  tabPanel("Miscellaneous")
                   )
  server <- function(input, output,session) {
    
    
    
    output$info <- renderPrint({
      nearPoints(final,input$plot_click,threshold = 10, maxpoints = 1,addDist = TRUE)
    })
    
    observe({
      updateSelectInput(session, "GENE", choices = as.character(final[final$SAMPLE==input$SAMPLE, "GENE"]%>% discard(is.na)))
    })
    
    observe({
      updateSelectInput(session, "SGDID", choices  = as.character(final[final$GENE==input$GENE, "SGDID"] %>% discard(is.na) %>% unique()))
    })
    
    output$url <- renderUI({
      url <- a("Learn about Gene",href=paste0(link,input$SGDID),class="btn btn-default")
      url
    })
    
    output$pdf_viewer <- renderUI({ tags$iframe(
      style="height:1000px;width:100%;scrolling=yes",
      src = "Black_box.pdf") }) 
    
    output$contents <- renderTable({
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "txt", "Please upload a VCF file"))
      
      final <- read.table(file$datapath, header = input$header)

      final %>% filter(SAMPLE==input$SAMPLE[1])
    })
    

    
    

    
    output$plot1 <- renderPlot({
      #df()
      final %>% 
        mutate(Chromosome=forcats::fct_relevel(Chromosome,'I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','M')) %>%
        ggplot() +
        facet_grid(vars(Chromosome),switch = "y") +
        geom_segment(aes(color=Chromosome),x = 1, y = 0, xend = final$Length, yend = 0, size=4.1,lineend = "round") +
        geom_segment(x = final$cent1, y = 0, xend = final$cent2, yend = 0, size=4.1,lineend = "round", color="black") +
        scale_color_manual(values=pl_palette("lorax",17)) +
        geom_point(aes(x=POS,y=0),shape = "|", size=2.9, data=final
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
      
      final %>% filter(SAMPLE==input$SAMPLE[1]) %>% 
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
      
      
      final %>% mutate(transition=paste(REF,"_",ALT, sep=""))  %>% 
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
    
    output$plot3 <- renderPlot({
      #df2()
      
      final %>% 
        filter(SAMPLE==input$SAMPLE[1]) %>%
        filter(GENE==input$GENE[1]) %>%
        mutate(ANNOTATION= gsub("'","",ANNOTATION)) %>%
        mutate(AA_POS = if_else(ANNOTATION=="5-upstream",-20,as.numeric(AA_POS))) %>%
        ggplot(aes(x=as.numeric(AA_POS),y=1)) + 
        #facet_grid(rows=vars(GENE))+
        geom_hline(yintercept=0, linetype=2,alpha=.2)+
        #geom_segment(aes(x=-10,xend=PROTEIN_LENGTH+10,y=0,yend=0), size=20, color = "pink") + 
        geom_segment(aes(x=0,xend=PROTEIN_LENGTH,y=0,yend=0), size=15, color = "cornflowerblue") +
        geom_segment(aes(x=as.numeric(AA_POS),xend=as.numeric(AA_POS),y=0,yend=1), color = "pink") +
        geom_point(aes(x=as.numeric(AA_POS), color=ANNOTATION),y=1, size=2)+
        ylim(c(-2,2))+ 
        geom_label_repel(aes(label = PROTEIN),
                         box.padding   = 1, 
                         point.padding = 1,
                         segment.color = 'grey50') +
        ggtitle(as.character(input$GENE[1]))+
        theme_classic() +
        theme(axis.title.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title = element_text(hjust = 0.5),
              axis.text.y = element_blank()) + 
        xlab("Amino acid position")
    })
    
    output$text <- renderText({ paste("Mutations Types:"," - A nonsynonymous substitution is a nucleotide mutation that alters the amino acid sequence of a protein.",
                                 "- A synonymous mutation is a change in the DNA sequence that codes for amino acids in a protein sequence," ,"but does not change the encoded amino acid.",
                                 "- The 5′ untranslated region (also known as 5′ UTR) is the region of a messenger RNA (mRNA) that is directly","upstream from the initiation codon. It is not usually translated.",
                                 "- Intergenic regions are the stretches of DNA located between genes.",
                                 "- An autonomously replicating sequence (ARS) contains the origin of replication in the yeast genome.", sep="\n")})
  }
  
  shinyApp(ui, server)

