# List of required packages
required_packages <- c("devtools", "devtools", "shinythemes", "DBI", "RSQLite", 
                       "ggplot2","dplyr", "tidyr", "forcats", "ggrepel",
                       "purrr", "PLColors", "shinyjs", "viridis", "plotly")

# Check if packages are installed, and if not, install them
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

#loading necessary libraries
library(devtools)
library(shiny)
library(shinythemes)
library(DBI)
library(RSQLite)
library(ggplot2) ## visualization data
library(dplyr)
library(viridis)
library(tidyr) ## handling data, df must be tidy to use a lot of packages 
library(forcats)
library(ggrepel)
library(purrr)
library(shinyjs)
library(plotly)



#loading in the VCF file to display initial choices, later turns into reactive val called mutation_data that includes manually updated data
mut_backend <- read.csv("all_yEvo_vcf.csv") #used to be called final, and used to run off of final_allVCF

#loading in the genes data file
genes_info <- read.csv("gene_info.csv")

#loading in the chromosomes data file
chrom_info <- read.csv("chromosome_info.csv")

#need to add this to upload the yEvo icon the theme 
addResourcePath(prefix = 'img', directoryPath = 'img')

#create a variable called link that stores the base SGD database for locus 
link = "https://www.yeastgenome.org/locus/"

#how will our layout look like for the app 
ui <-  navbarPage(
  useShinyjs(),
  title = div(img(src="img/yEvo_logo.png",
                  filetype = "image/png",
                  style="margin-top: -14px;
                         padding-right:10px;
                         padding-bottom:10px",
                  height = 60)
  ),
  
  windowTitle="yEvo",
  theme = shinytheme("cerulean"),
  tabPanel("Data Visualizations",
           sidebarLayout(
             # Left side, Class vs Cumulative View and options 
             sidebarPanel(
               fileInput("datafile", "Optional: Upload additional CSV File", accept = ".csv"),
               # conditionalPanel(
               #   # Asks for instructor/year so we can modify user uploaded file so we can combine it
               # with our current data
               #   condition = "output.filesUploaded",
               #   textInput("inputted_instructor", "Who is your Instructor"),
               #   textInput("inputted_year", "What is the current year")
               # ),
               radioButtons("View", "Select an option:",
                            choices = c("View By Class", "View By Selection Condition"),
                            selected = character(0)),
               div("", style = "height: 10px;"),  # Create a 10px vertical space
               
               #create download button here
               downloadButton("downloadBtn", "Download the following data"),
               
               div("", style = "height: 10px;"),  # Create a 10px vertical space
               
               # Only shows on condition observeEvent
               conditionalPanel(
                 # links condition to button via button key 
                 condition = "input.uploadData || input.classView || output.selectedClassView",
                 selectInput("instructor", "Instructor", choices = c('')),
                 selectInput("year", "Year", choices = c('')),
                 # PLEASE NOTE: "sample" is called "Sample Name" within ui to make it easier to understand, but the official name in the df is sample                 
                 selectInput("sample", "Sample Name", choices = c('')),
               ),
               conditionalPanel(
                 condition = "input.cumulView || output.selectedCumulView",
                 selectInput("condition", "Condition", choices = c('All Selected', mut_backend %>% count(condition) %>% pull(condition))),
                 selectInput("background", "Background", choices = c('')),
               ),
             ),
             # Right side, Data Visualization
             mainPanel(
               tabsetPanel(
                 type = "tabs",
                 #tabPanel("Chromosome Map", plotlyOutput("chromPlot", brush = brushOpts(id = "plot_brush", fill = "#ccc", direction = "x")),verbatimTextOutput("info")),
                 tabPanel("Chromosome Map", plotlyOutput("chromPlot"),verbatimTextOutput("info")),
                 tabPanel("Variant Pie Chart", plotOutput("varPieChart", click = "plot_click"), verbatimTextOutput("text")),
                 tabPanel("SNP Counts", plotOutput("snpCountPlot", click = "plot_click")),
                 tabPanel("Gene View", value = "Geneview", plotOutput("geneViewPlot", dblclick = "geneViewPlot_dblclick", brush = brushOpts(id = "geneViewPlot_brush", resetOnNew = TRUE)),
                          selectInput("GENE", "Gene", choices = c('')),
                          uiOutput("url"),
                          verbatimTextOutput("text1")),
                 tabPanel("Table", tableOutput("data_table")),
               )
             )
           )
  ),
  tabPanel("Tutorial",
           div(img(src="img/how-to-WIDE.png",
                   height="65%", width="65%"),
               style="text-align: center;")
  ),
  tabPanel("Background",
           uiOutput("pdf_viewer"))
) #END OF UI

# Now entering server, which handles everything dynamically
server <- function(input, output,session) {
  #initially setting default file of all mutation data
  mutation_data <- reactiveVal(read.csv("all_yEvo_vcf.csv")) 
  shinyjs::hide("cumulDropdowns") # Initially hide cumulative drop downs
  
  # Displays Chromosome Map info; filtering by sample
  output$info <- renderText({
    return("Drag over variant tick mark to see details\n")
  })
  # Initialize a reactive variable for the dataframe
  # Function to read and append the uploaded data to the cumulative dataframe
  #TODO:change "final" name to new file name (new file as in the new system we are making)
  observeEvent(input$datafile, {
    file <- input$datafile
    if (!is.null(file) && all(names(mut_backend) == names(read.csv(file$datapath, sep = ",")))) {
      df <- rbind(mut_backend,read.csv(file$datapath, sep = ","))
      #df <- read.csv(file$datapath, sep = ",")
      mutation_data(df)
    } else {
      mutation_data(mut_backend)
    }
  })
  
  
  #filtering dataframe based on menu selection
  filtered_data <- reactive({
    data <- mutation_data()
    if(!is.null(input$View)){
      if (input$View == "View By Class") {
        # Get the selected values from the dropdown menus
        selected_instructor <- input$instructor
        selected_year <- input$year
        selected_sample <- input$sample
        #filtering based on selections if NOT all selected
        if (selected_instructor != "All Selected") {
          data <- data %>% filter(instructor == selected_instructor)
        }
        if (selected_year != "All Selected") {
          data <- data %>% filter(year == selected_year)
        }
        if (selected_sample != "All Selected") {
          data <- data %>% filter(sample == selected_sample)
        }
      } else if (input$View == "View By Selection Condition") {
        selected_condition <- input$condition
        selected_background <- input$background
        
        if (selected_condition != "All Selected") {
          data <- data %>% filter(condition == selected_condition)
        }
        
        if (selected_background != "All Selected") {
          data <- data %>% filter(background == selected_background)
        }
      }
    }
    data 
  })
  
  #download button functionality
  output$downloadBtn <- downloadHandler(
    filename = function() {
      # Set the filename for the downloaded file
      if(is.null(input$View)){
        "master_table.csv"
      }else if (input$View == "View By Class") {
        selected_instructor <- input$instructor
        selected_year <- input$year
        selected_sample <- input$sample
        if (selected_sample != "All Selected") {
          paste0(selected_instructor, "_",selected_year,"_",selected_sample,".csv")
        }else if (selected_year != "All Selected") {
          paste0(selected_instructor, "_",selected_year,".csv")
        }else if (selected_instructor != "All Selected") {
          paste0(selected_instructor,".csv")
        }
        
      }else if (input$View == "View By Selection Condition") {
        if(input$background != "All Selected"){
          paste0(input$condition,"_",input$background,".csv")
        }else{
          paste0(input$condition,".csv")
        }
      }
    },
    content = function(file) {
      # Write the data to a CSV file
      write.csv(filtered_data(), file)
    }
  )
  
  # storing a bool to see if a file has been uploaded
  # if a file has be uploaded, using the condition that if output$filesUploaded 
  # is true, we can auto-open the class view
  output$filesUploaded <- reactive({
    val <- !(is.null(input$datafile))
  })
  
  outputOptions(output, 'filesUploaded', suspendWhenHidden=FALSE)
  
  output$selectedClassView <- reactive({
    if(!is.null(input$View)){
      value <- (input$View == "View By Class")
    }
  })
  outputOptions(output, 'selectedClassView', suspendWhenHidden=FALSE)
  
  output$selectedCumulView <- reactive({
    if(!is.null(input$View)){
      value <- (input$View == "View By Selection Condition")
    }
  })
  outputOptions(output, 'selectedCumulView', suspendWhenHidden=FALSE)
  
  # Render the dataframe in the tableOutput
  output$data_table <- renderTable({
    filtered_data() %>%
      select(CHROM, POS, ANNOTATION, GENE, PROTEIN, condition, instructor, year, sample, REF, ALT)
  })
  
  #Display settings
  observe({
    if(!is.null(input$View)){
      if (input$View == "View By Class") { 
        shinyjs::disable("cumulView")
        shinyjs::enable("classView")
      }
      else if(input$View == "View By Selection Condition"){
        shinyjs::enable("cumulView")
        shinyjs::disable("classView")
        #shinyjs:disable("background")
      }
    } else {
      shinyjs::disable("classView")
      shinyjs::disable("cumulView")
    }
  })
  
  #Handling behaviors for button selections
  observe({
    options <- c(as.character(mutation_data() %>% filter(condition == input$condition) %>% pull(background)))
    print(options)
    if(length(unique(options)) != 1){ 
      updateSelectInput(session, "background", choices = c('All Selected', options))
      shinyjs::enable("background")
    }else{
      updateSelectInput(session, "background", choices = c(options))
      shinyjs::disable("background")
    }
  })
  #making the background not a button before selecting a condition
  observe({
    if(input$condition == "All Selected"){
      shinyjs::disable("background")
    }
  })
  
  observe({
    updateSelectInput(session, "instructor", choices = c("All Selected", unique(mutation_data()$instructor)))
  })
  
  
  observe({
    if (input$instructor == "All Selected") {
      updateSelectInput(session, "year", choices = c("All Selected", unique(mutation_data()$year)))
    } else {
      updateSelectInput(session, "year", choices = c("All Selected", as.character(mutation_data()[mutation_data()$instructor == input$instructor, "year"])))
    }
  })
  
  observe({
    updateSelectInput(session, "instructor", choices = c('All Selected', unique(mutation_data()$instructor)))
  }) 
  
  
  observe({
    updateSelectInput(session, "sample", choices = c("All Selected", as.character(mutation_data() %>% filter(instructor==input$instructor) %>% filter(year==input$year) %>% pull(sample))))
  })
  
  
  observe({
    updateSelectInput(session, "GENE", choices = if(input$sample!="All Selected") { as.character(mutation_data()[mutation_data()$sample==input$sample, "GENE"]%>% discard(is.na))
    } else {as.character(filtered_data() %>% pull(GENE) %>% discard(is.na)) })
    
  })
  
  #TODO: make sure gene info and sgdid info is coming from the right place ->should all be handled within filtering function(?)
  # Learn about Gene button within gene viewer
  sgdid <- reactiveValues(value = NULL)
  output$url <- renderUI({
    sgdid_values <- as.character(genes_info[genes_info$GENE == input$GENE, "SGDID"] %>% discard(is.na) %>% unique())
    sgdid$value <- sgdid_values
    url <- a("Learn about Gene",href=paste0(link, sgdid$value),class="btn btn-default", target='_blank')
    url
  })
  
  output$pdf_viewer <- renderUI({ tags$iframe(
    style="height:1000px;width:100%;scrolling=yes",
    src = "Black_box.pdf") }) 
  
  tabPanel("Background",
           uiOutput("pdf_viewer") )
  
  #TODO: replace chromosome and gene data w the correct path
  
  # Define the desired order of categories
  desired_order <- c('chrM', 'chrXVI', 'chrXV', 'chrXIV', 'chrXIII', 'chrXII', 'chrXI', 'chrX', 'chrIX', 'chrVIII', 'chrVII', 'chrVI', 'chrV', 'chrIV', 'chrIII', 'chrII','chrI')
  # Convert category to a factor with the desired order
  chrom_info$CHROM <- factor(chrom_info$CHROM, levels = desired_order)

  #merge(genes_info,chromosomes, by = “CHROM”)
  # Define a mapping from chromosome names to numbers
  chromosome_mapping <- c(
    'chrM' = 1,
    'chrXVI' = 2,
    'chrXV' = 3,
    'chrXIV' = 4,
    'chrXIII' = 5,
    'chrXII' = 6,
    'chrXI' = 7,
    'chrX' = 8,
    'chrIX' = 9,
    'chrVIII' = 10,
    'chrVII' = 11,
    'chrVI' = 12,
    'chrV' = 13,
    'chrIV' = 14,
    'chrIII' = 15,
    'chrII' = 16,
    'chrI' = 17
  )
  mutations_locations <- reactive({
    # Extract the current values of the reactive data frames
    mutation_data_value <- filtered_data()
    # Merge the data frames based on the “REGION” column
    merged <- merge(mutation_data_value, genes_info, by = 'REGION')
    # Adding chrom as num
    merged$chromosome_as_num <-chromosome_mapping[merged$CHROM.x]
    #returning merged
    merged
  })
  
  # Chromosomal data for chromoplot
  chromosomes <- data.frame(
    chromosome = chrom_info$CHROM,
    length = chrom_info$length
  )
  
  #To pull out only the genes that mutated
  mutated_genes <- reactive({
    # To pull out only the genes that mutated
    data.frame(
      chromosome = mutations_locations()$CHROM.x,
      start = mutations_locations()$START,
      end = mutations_locations()$STOP,
      geneName = mutations_locations()$GENE.y,
      chrom_as_num = mutations_locations()$chromosome_as_num
    )
  })
  
  # Create an empty dataframe to store the final results
  final_gene <- reactive({
  final_gene_static <- data.frame()
    
  # Iterate through unique genes
  unique_genes <- unique(filtered_genes_data$GENE.x)
  for (gene in unique_genes) {
    # Filter data for the current gene
    gene_data <- filtered_genes_data %>% filter(GENE.x == gene)
    
    # Count the number of repetitions
    num_repeats <- nrow(gene_data)
    
    # Get the associated columns
    gene_info <- gene_data[1, c("CHROM.x", "START", "STOP", "GENE.x")]
    
    # Add the number of repeats as a new column
    gene_info$repeats <- num_repeats
    
    # Add this gene to the final dataframe
    final_gene_static <- rbind(final_gene_static, gene_info)
  }
  final_gene_static$chrom_as_num <- chromosome_mapping[final_gene_static$CHROM.x]
  
  final_gene_static
  })
  
  
  #chromosomes$chromosome_as_num <- chromosome_mapping[chromosomes$chromosome]
  output$chromPlot <- renderPlotly({
    print(mutated_genes())
    print(chromosomes)
    # Plotting
    p <- ggplot() +
      geom_bar(data = chromosomes, aes(x = length, y = chromosome), stat = 'identity', fill = 'lightblue', width = 0.5) + # swapped x and y
      geom_rect(data = final_gene(), aes(ymin = chrom_as_num - 0.4, # swapped ymin and ymax
                                            ymax = chrom_as_num + 0.4,
                                            xmin = START,
                                            xmax = STOP,
                text = paste("Gene Name: ",GENE.x),
                fill = repeats), 
    
    , alpha = 0.5) +
    scale_fill_gradient(low = "navy", high = "red",limits = c(0, 25)) +
      labs(title = 'Location of mutations along chromosomes',
           y = 'Chromosome', # changed x-axis label to Chromosome
           x = 'Position along chromosome') # changed y-axis label to Length
    # Convert ggplot2 plot to plotly
    p <- ggplotly(p)
    # Add formatting
    p <- layout(
      p,
      plot_bgcolor = 'rgba(0,0,0,0)',   # Set plot background color to transparent
      paper_bgcolor = 'rgba(0,0,0,0)',  # Set paper background color to transparent
      xaxis = list(showgrid = FALSE),  # Remove x-axis gridlines
      yaxis = list(showgrid = FALSE)   # Remove y-axis gridlines
    )
  })
  
  
  
  # Maps specific annotation to specific color
  annotat_colormap <- c()
  
  output$varPieChart <- renderPlot({
    # Color vector for each annotation in Pie Chart
    # just add the same number of colors as number of annotations
    # ex. if there are 10 unique annotations, put 10 unique colors
    # in this color vector
    color_vector <- c("#9edae5", "#17becf", "#dbdb8d", "#bcbd22", "#c7c7c7",
                      "#7f7f7f", "#f7b6d2", "#e377c2", "#c49c94", "#8c564b",
                      "#c5b0d5", "#9467bd", "#ff9896", "#d62728", "#98df8a",
                      "#2ca02c", "#ffbb78", "#ff7f0e", "#aec7e8", "#1f77b4")
    
    # gives us the number of unique annotations in filtered data
    unique_annotations <- filtered_data() %>%
      distinct(ANNOTATION) %>% pull(ANNOTATION)
    
    # if there are new annotations that have not been mapped to a color
    # put them in this named vector with [annotation = color] 
    new_anno_vector <- c()
    # keep track of which unique index in the color_vector we will choose
    # for our new annotations
    cur_anno_length = length(annotat_colormap) + 1
    # check along all unique annotations
    for (i in seq_along(unique_annotations)) {
      # if there are new annotations, add it to the color map 
      # and set it to next available color
      if (!(unique_annotations[i] %in% names(annotat_colormap))){
        # get the new unused color for our new annotation
        new_color = color_vector[cur_anno_length]
        # add this new annotation and color name-value pair to a vector
        new_anno_vector <- c(new_anno_vector, setNames(new_color,unique_annotations[i]))
        # increment to cur_anno_length to next unique color
        cur_anno_length = cur_anno_length + 1
      }
    }
    # now combine back into annotat_colormap
    # and permanently update (does not reset unless you close and rerun the program)
    annotat_colormap <<- c(annotat_colormap, new_anno_vector)
    
    blank_theme <- theme_minimal()+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold"))
    
    filtered_data() %>% 
      count(ANNOTATION) %>% 
      mutate(percent=n/sum(n)*100) %>% 
      ggplot(aes(x="",y=percent,fill=ANNOTATION)) + 
      geom_bar(stat="identity", width=1) +
      coord_polar("y", start=0) + 
      ggtitle("Percentage of Variants by Type") + 
      geom_text(aes(label = round(percent), digits = 0),position = position_stack(vjust = 0.5),col="white") + 
      blank_theme + 
      scale_fill_manual(values=annotat_colormap)
  })
  
  
  #TODO: figure out what ref, alt, are from and what info we need here?
  output$snpCountPlot <- renderPlot({
    # counting number of transitions to display to set for colors and plot
    num <- filtered_data() %>% mutate(transition=paste(REF,"_",ALT, sep="")) %>% 
      mutate(length = nchar(transition)) %>% mutate(transition = if_else(nchar(transition) > 3,"Indel",transition)) %>% 
      count(transition) %>% summarise(n = n()) %>% as.numeric()
    
    filtered_data() %>% mutate(transition=paste(REF,"_",ALT, sep=""))  %>% 
      mutate(length = nchar(transition)) %>% 
      mutate(transition = if_else(nchar(transition) > 3,"Indel",transition)) %>%
      ggplot(aes(x=as.factor(transition),fill=as.factor(transition))) + 
      geom_bar() + theme_bw() + 
      scale_fill_manual(values=viridis(n = num, begin = 0.4, end = 1)) + 
      theme(legend.position = "none",
            strip.text.y.left = element_text(angle = 0),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=.5)) + 
      ggtitle("Single Nucleotide transitions")  + xlab("SNP call")
  })
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  #TODO: AA_POS what is going on here lol and do we have the right info
  output$geneViewPlot <- renderPlot({
    xlength <- genes_info %>%
      filter(GENE==input$GENE) %>% pull(PROTEIN_LENGTH) %>% unique() %>% as.numeric()
    
    filtered_data() %>% mutate("AA_POS" = stringr::str_extract(PROTEIN, "([0-9])+")) %>% 
      filter(GENE==input$GENE) %>%
      mutate(ANNOTATION= gsub("'","",ANNOTATION)) %>%
      mutate(AA_POS = if_else(ANNOTATION=="5-upstream",-15,as.numeric(AA_POS))) %>%
      ggplot(aes(x=as.numeric(AA_POS),y=.5)) + 
      geom_hline(yintercept=0, linetype=2,alpha=.2)+
      geom_segment(aes(x=0,xend=xlength,y=0,yend=0), size=15, color = "cornflowerblue") +
      geom_segment(aes(x=as.numeric(AA_POS),xend=as.numeric(AA_POS),y=0,yend=.5), color = "pink") +
      geom_point(aes(x=as.numeric(AA_POS), color=ANNOTATION),y=0.5, size=2)+
      ylim(c(-0.2, 1.2))+ 
      xlim(-50,xlength)+
      geom_text_repel(aes(label = PROTEIN),
                      box.padding   = 2, 
                      point.padding = 1,
                      segment.color = 'grey50',
                      min.segment.length = 0
      ) +
      ggtitle(as.character(input$GENE))+
      theme_classic(base_size=18) +
      theme(axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_blank()) + 
      xlab("Amino acid position") +
      theme(axis.text.x = element_text(size = 8)) +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) + 
      theme(plot.margin = margin(0, 0, 0, 0)) + 
      annotate(
        "text", x = 1, y = Inf, label = "Drag over mutations to see more",
        hjust = 0, vjust = 2, color = "black", size = 5
      )
    
  }, width = 750)
  
  #TODO: how can we find the y values to properly match up with each chromosome?
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$geneViewPlot_dblclick, {
    brush <- input$geneViewPlot_brush
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
  
  observeEvent(input$append_btn, {
    new_csv_path <- input$new_csv$datapath
    
  }) 
}

shinyApp(ui, server)