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

#final <- read.table("final_MASTERVCF.txt", header=TRUE)
final <- read.csv("final_allVCF.csv")

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
        #fileInput("file1", "Choose VCF File", accept = ".txt"),
        #checkboxInput("header", "Header", TRUE),
        selectInput("instructor","Instructor",
                    choices = final %>% count(instructor) %>% pull(instructor)),
        
        selectInput("year","Year",
                    choices = c('')),
        
        selectInput("sample", "Sample",
                    choices = c('')),
        #br(),
        
        h5("If Sample = None Selected, you can browse visualizations for all data based on condition and background combination chosen "),
        
        selectInput("condition","Condition",
                    choices = final %>% count(condition) %>% pull(condition)),
        
        selectInput("background","Background",
                    choices = c('')),
        
        h5("Gene and SGDID menu will update either based on the sample chosen or the condition/background combination chosen."),
        
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
                    tabPanel("Gene View", value="Geneview",plotOutput("plot3", dblclick = "plot3_dblclick",
                                                                      brush = brushOpts(
                                                                        id = "plot3_brush",
                                                                        resetOnNew = TRUE)),verbatimTextOutput("text1")),
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
    
    
    
   # output$info <- renderPrint({
  #    nearPoints(final,input$plot_click,threshold = 10, maxpoints = 1,addDist = TRUE)
  #  })
    
    
    output$info <- renderText({
      xy_range_str <- function(e) {
        if(is.null(e)) return("Drag over variant tick mark to see details\n")
        paste0("Variant Gene: ",final %>% filter(if (input$sample != "None Selected") {sample==input$sample} else {condition==input$condition}) %>% filter(POS > round(e$xmin, 1)) %>% filter(POS < round(e$xmax, 1)) %>% select(GENE), "\n",
              "Reference: ", final %>% filter(if (input$sample != "None Selected") {sample==input$sample} else {condition==input$condition}) %>% filter(POS > round(e$xmin, 1)) %>% filter(POS < round(e$xmax, 1)) %>% select(REF),"\n",
              "Variant: ",final %>% filter(if (input$sample != "None Selected") {sample==input$sample} else {condition==input$condition}) %>% filter(POS > round(e$xmin, 1)) %>% filter(POS < round(e$xmax, 1)) %>% select(ALT),"\n",
              "Position: ",final %>% filter(if (input$sample != "None Selected") {sample==input$sample} else {condition==input$condition}) %>% filter(POS > round(e$xmin, 1)) %>% filter(POS < round(e$xmax, 1)) %>% select(POS),"\n",
              "Chromosome: ",final %>% filter(if (input$sample != "None Selected") {sample==input$sample} else {condition==input$condition}) %>% filter(POS > round(e$xmin, 1)) %>% filter(POS < round(e$xmax, 1)) %>% select(Chromosome))
      }
      
      paste0(
        xy_range_str(input$plot_brush)
      )
    })
    
    
    observe({
      updateSelectInput(session, "year", choices = as.character(final[final$instructor==input$instructor, "year"]))
    })
    
    observe({
      updateSelectInput(session, "sample", choices = c("None Selected", as.character(final %>% filter(instructor==input$instructor) %>% filter(year==input$year) %>% pull(sample))))
    })
    
    observe({
      updateSelectInput(session, "background", choices = as.character(final %>% filter(condition==input$condition) %>% pull(background)))
    })
    
    observe({
      updateSelectInput(session, "GENE", choices = if(input$sample!="None Selected") { as.character(final[final$sample==input$sample, "GENE"]%>% discard(is.na))
        } else {as.character(final[final$condition==input$condition, "GENE"]%>% discard(is.na)) })
    })
    
    observe({
      updateSelectInput(session, "SGDID", choices  = as.character(final[final$GENE==input$GENE, "SGDID"] %>% discard(is.na) %>% unique()))
    })
    
    output$url <- renderUI({
      url <- a("Learn about Gene",href=paste0(link,input$SGDID),class="btn btn-default", target='_blank')
      url
    })
    
    output$pdf_viewer <- renderUI({ tags$iframe(
      style="height:1000px;width:100%;scrolling=yes",
      src = "Black_box.pdf") }) 
    
    output$contents <- renderTable({
      #file <- input$file1
      #ext <- tools::file_ext(file$datapath)
      
      #req(file)
      #validate(need(ext == "csv", "Please upload a VCF file"))
      
      #final <- read.csv(file$datapath, header = input$header)
      if(input$sample!="None Selected") {
        final %>% filter(sample==input$sample)
      }else {
        final %>% filter(condition==input$condition)
      }
      
    })
    
    
    output$plot1 <- renderPlot({
      if (input$sample != "None Selected") {
      final %>% 
        mutate(Chromosome=forcats::fct_relevel(Chromosome,'I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','M')) %>%
        ggplot() +
        facet_grid(vars(Chromosome),switch = "y") +
        geom_segment(aes(color=Chromosome),x = 1, y = 0, xend = final$Length, yend = 0, size=4.1,lineend = "round") +
        geom_segment(x = final$cent1, y = 0, xend = final$cent2, yend = 0, size=4.1,lineend = "round", color="black") +
        scale_color_manual(values=pl_palette("lorax",17)) +
        geom_point(aes(x=POS,y=0),shape = "|", size=2.9, data=final
                   %>% mutate(Chromosome=forcats::fct_relevel(Chromosome,'I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','M'))%>% 
                     filter(sample==input$sample)) + 
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
      } else {
        final %>% 
          mutate(Chromosome=forcats::fct_relevel(Chromosome,'I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','M')) %>%
          ggplot() +
          facet_grid(vars(Chromosome),switch = "y") +
          geom_segment(aes(color=Chromosome),x = 1, y = 0, xend = final$Length, yend = 0, size=4.1,lineend = "round") +
          geom_segment(x = final$cent1, y = 0, xend = final$cent2, yend = 0, size=4.1,lineend = "round", color="black") +
          scale_color_manual(values=pl_palette("lorax",17)) +
          geom_point(aes(x=POS,y=0),shape = "|", size=2.9, data=final
                     %>% mutate(Chromosome=forcats::fct_relevel(Chromosome,'I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','M'))%>% 
                       filter(condition==input$condition) %>% filter(background==input$background)) + 
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
      }
      
      })
    
    output$plot <- renderPlot({
      if (input$sample != "None Selected") {
        num <- final %>% filter(condition==input$condition) %>% 
          filter(background==input$background) %>%
          count(ANNOTATION) %>% summarise(n = n()) %>% as.numeric()
        
      blank_theme <- theme_minimal()+
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title=element_text(size=14, face="bold"))
      
      final %>% filter(sample==input$sample) %>% 
        count(ANNOTATION) %>% 
        mutate(percent=n/sum(n)*100) %>% 
        ggplot(aes(x="",y=percent,fill=ANNOTATION)) + 
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) + 
        ggtitle("Percentage of Variants by Type") + 
        geom_text(aes(label = round(percent), digits = 0),position = position_stack(vjust = 0.5),col="white") + 
        blank_theme + 
        scale_fill_manual(values=pl_palette("lorax",num))
      } else {
         num <- final %>% filter(condition==input$condition) %>% 
          filter(background==input$background) %>%
          count(ANNOTATION) %>% summarise(n = n()) %>% as.numeric()
           
        blank_theme <- theme_minimal()+
          theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title=element_text(size=14, face="bold"))
        
        final %>% filter(condition==input$condition) %>% 
          filter(background==input$background) %>%
          count(ANNOTATION) %>% 
          mutate(percent=n/sum(n)*100) %>% 
          ggplot(aes(x="",y=percent,fill=ANNOTATION)) + 
          geom_bar(stat="identity", width=1) +
          coord_polar("y", start=0) + 
          ggtitle("Percentage of Variants by Type") + 
          geom_text(aes(label = round(percent), digits = 0),position = position_stack(vjust = 0.5),col="white") + 
          blank_theme + 
          scale_fill_manual(values=pl_palette("lorax",num))
      }
    })
    
    
    output$plot2 <- renderPlot({
      if(input$sample!="None Selected") {
      num <- final %>% mutate(transition=paste(REF,"_",ALT, sep=""))  %>% select(transition,sample) %>% mutate(length = nchar(transition)) %>% 
        filter(sample==input$sample[1]) %>%
        count(transition) %>% 
        summarise(n = n()) %>% as.numeric()
      
      
      final %>% mutate(transition=paste(REF,"_",ALT, sep=""))  %>% 
        select(transition,sample) %>% 
        filter(sample==input$sample[1]) %>%
        mutate(length = nchar(transition)) %>% 
        #filter(length >= 3) %>%  
        mutate(transition = if_else(nchar(transition) > 3,"Indel",transition)) %>%
        ggplot(aes(x=as.factor(transition),fill=as.factor(transition))) + 
        geom_bar() + theme_bw() + 
        scale_fill_manual(values=pl_palette("lorax",num)) + 
        theme(legend.position = "none",
              strip.text.y.left = element_text(angle = 0),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=.5)) + 
        ggtitle("Single Nucleotide transitions")  + xlab("SNP call")
      } else {
        num <- final %>% mutate(transition=paste(REF,"_",ALT, sep=""))  %>% select(transition,sample,condition,background) %>% mutate(length = nchar(transition)) %>% 
          filter(condition==input$condition) %>%
          filter(background==input$background) %>%
          count(transition) %>% 
          summarise(n = n()) %>% as.numeric()
        
        
        final %>% mutate(transition=paste(REF,"_",ALT, sep=""))  %>% 
          select(transition,sample,condition,background) %>% 
          filter(condition==input$condition) %>%
          filter(background==input$background) %>%
          mutate(length = nchar(transition)) %>% 
          #filter(length >= 3) %>%  
          mutate(transition = if_else(nchar(transition) > 3,"Indel",transition)) %>%
          ggplot(aes(x=as.factor(transition),fill=as.factor(transition))) + 
          geom_bar() + theme_bw() + 
          scale_fill_manual(values=pl_palette("lorax",num)) + 
          theme(legend.position = "none",
                strip.text.y.left = element_text(angle = 0),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=.5)) + 
          ggtitle("Single Nucleotide transitions")  + xlab("SNP call")
      }
    })
    
    ranges <- reactiveValues(x = NULL, y = NULL)
    
    output$plot3 <- renderPlot({
    if(input$sample != "None Selected") {
      
      xlength <- final %>% filter(condition==input$condition) %>%
        filter(background==input$background) %>%
        filter(GENE==input$GENE[1]) %>% pull(PROTEIN_LENGTH) %>% unique() %>% as.numeric()
      
      final %>% 
        filter(sample==input$sample[1]) %>%
        filter(GENE==input$GENE[1]) %>%
        mutate(ANNOTATION= gsub("'","",ANNOTATION)) %>%
        mutate(AA_POS = if_else(ANNOTATION=="5-upstream",-15,as.numeric(AA_POS))) %>%
        ggplot(aes(x=as.numeric(AA_POS),y=1)) + 
        #facet_grid(rows=vars(GENE))+
        geom_hline(yintercept=0, linetype=2,alpha=.2)+
        #geom_segment(aes(x=-10,xend=PROTEIN_LENGTH+10,y=0,yend=0), size=20, color = "pink") + 
        geom_segment(aes(x=0,xend=PROTEIN_LENGTH,y=0,yend=0), size=15, color = "cornflowerblue") +
        geom_segment(aes(x=as.numeric(AA_POS),xend=as.numeric(AA_POS),y=0,yend=1), color = "pink") +
        geom_point(aes(x=as.numeric(AA_POS), color=ANNOTATION),y=1, size=2)+
        ylim(c(-2,2))+ 
        xlim(-20,xlength)+
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
        xlab("Amino acid position") +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    } else {
      
      xlength <- final %>% filter(condition==input$condition) %>%
        filter(background==input$background) %>%
        filter(GENE==input$GENE[1]) %>% pull(PROTEIN_LENGTH) %>% unique() %>% as.numeric()
      
      final %>% 
        filter(condition==input$condition) %>%
        filter(background==input$background) %>%
        filter(GENE==input$GENE[1]) %>%
        mutate(ANNOTATION= gsub("'","",ANNOTATION)) %>%
        mutate(AA_POS = if_else(ANNOTATION=="5-upstream",-15,as.numeric(AA_POS))) %>%
        ggplot(aes(x=as.numeric(AA_POS),y=1)) + 
        #facet_grid(rows=vars(GENE))+
        geom_hline(yintercept=0, linetype=2,alpha=.2)+
        #geom_segment(aes(x=-10,xend=PROTEIN_LENGTH+10,y=0,yend=0), size=20, color = "pink") + 
        geom_segment(aes(x=0,xend=PROTEIN_LENGTH,y=0,yend=0), size=15, color = "cornflowerblue") +
        geom_segment(aes(x=as.numeric(AA_POS),xend=as.numeric(AA_POS),y=0,yend=1), color = "pink") +
        geom_point(aes(x=as.numeric(AA_POS), color=ANNOTATION),y=1, size=2)+
        ylim(c(-2,2))+ 
        xlim(-20,xlength)+
        geom_label_repel(aes(label = as.character(PROTEIN)),
                        # box.padding   = 1, 
                         #point.padding = 1,
                         segment.color = 'grey50') +
        ggtitle(as.character(input$GENE[1]))+
        theme_classic() +
        theme(axis.title.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title = element_text(hjust = 0.5),
              axis.text.y = element_blank()) + 
        xlab("Amino acid position") +
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    }
      
    })
    
    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observeEvent(input$plot3_dblclick, {
      brush <- input$plot3_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    
    output$text <- renderText({ paste("Mutations Types:","- A nonsynonymous substitution is a nucleotide mutation that alters the amino acid sequence of a protein.",
                                 "- A synonymous mutation is a change in the DNA sequence that codes for amino acids in a protein sequence," ,"but does not change the encoded amino acid.",
                                 "- The 5′ untranslated region (also known as 5′ UTR) is the region of a messenger RNA (mRNA) that is directly","upstream from the initiation codon. It is not usually translated.",
                                 "- Intergenic regions are the stretches of DNA located between genes.",
                                 "- An autonomously replicating sequence (ARS) contains the origin of replication in the yeast genome.", sep="\n")})
    
    output$text1 <- renderText({ "The plot can be zoomed in by clicking and draggin and then double-clicking on the box.\n Reset view by double clicking again."})
  }
  
  shinyApp(ui, server)

