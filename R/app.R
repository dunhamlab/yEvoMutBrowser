# If you are changing the 'masterfile' that the app is built off of (aka the initial database),
# This needs to be changed in two locations:
# mut_backend before the server is called, and mutation_data right after the server is called.

### Imports
#loading necessary libraries
library(DBI)
library(devtools)
library(dplyr)
library(forcats)
library(ggplot2) ## visualization data
library(ggrepel)
library(plotly)
library(purrr)
library(RSQLite)
library(scales)
library(shiny)
library(shinyjs)
library(shinythemes)
library(stringr)
library(tidyr) ## handling data, df must be tidy to use a lot of packages 
library(viridis)

yEvoMutBrowser <- function(...) {
  # loading in the VCF file to display initial choices, later turns into reactive val called mutation_data that includes manually updated data
  # THIS SHOULD NOT BE CHANGED IN THE CODE. If the overall master shifts you can modify it here, 
  # but mut_backend should not be written to anywere in the code it is designed to just fill in in the beginning, 
  # right before the reactive frame gets created.
  mut_backend <- read.csv("all_yEvo_vcf.csv") 

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
              selectionPanel(mut_backend),
              # Right side, Data Visualization
              mainPanel(
                tabsetPanel(
                  type = "tabs",
                  chromMapUI("chromMap"),
                  variantsUI("variants"),
                  snpCountUI("snpCount"),
                  geneViewUI("geneView"),
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
    tabPanel("How Does Sequencing Work?",
            uiOutput("pdf_viewer"))
  ) #END OF UI

  # Now entering server, which handles everything dynamically
  server <- function(input, output,session) {
    
    #initially setting default file of all mutation data
    #CHANGE MASTERFILE HERE IF NEEDED
    mutation_data <- reactiveVal(read.csv("all_yEvo_vcf.csv")) 
    
    shinyjs::hide("cumulDropdowns") # Initially hide cumulative drop downs
    
    # Displays Chromosome Map info; filtering by sample
    output$info <- renderText({
      return("Drag over variant tick mark to see details\n")
    })
    # Initialize a reactive variable for the dataframe
    # Function to read and append the uploaded data to the cumulative dataframe

    observeEvent(input$submit_teach_year, {
      file <- input$datafile
      data <- read.csv(file$datapath)
      required_columns <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "ANNOTATION", "REGION", "GENE", "PROTEIN", "seq_file", 
                            "background", "condition", "sample")
      if (!is.null(file) && all(required_columns %in% colnames(data))) {
        # Add instructor and year columns
        data$instructor <- rep(input$inputted_instructor, nrow(data))
        data$year <- rep(input$inputted_year, nrow(data))
        # Rearrange columns
        data <- data[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "ANNOTATION", "REGION", "GENE", "PROTEIN", "seq_file", 
                        "background", "condition", "instructor", "year", "sample")]
        df <- rbind(mut_backend, data)
        # df <- read.csv(file$datapath, sep = ",")
        mutation_data(df)

      } else {
        # Handle the case where required columns are missing
        print("Some required columns are missing.")
        mutation_data(mut_backend)
      }
    })
    
    #filtering dataframe based on menu selection, most things from here on out should be based on filtered_data()
    filtered_data <- reactive({
      data <- mutation_data()
          #filtering based on selections if NOT all selected
          if (input$instructor != "All Selected") {
            data <- data %>% filter(instructor == input$instructor)
          }
          if (input$year != "All Selected") {
            data <- data %>% filter(year == input$year)
          }
          if (input$sample != "All Selected") {
            data <- data %>% filter(sample == input$sample)
          }
          selected_condition <- input$condition
          selected_background <- input$background
          
          if (input$condition != "All Selected") {
            data <- data %>% filter(condition == input$condition)
          }
          if (input$background != "All Selected") {
            data <- data %>% filter(background == input$background)
          }
      data 
    })
    
    #download button functionality
    output$downloadBtn <- downloadHandler(
      filename = function() {
        # Set the filename for the downloaded file
        if(is.null(input$View) || (input$instructor == "All Selected" && input$condition == "All Selected")){
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
          updateTextInput(session, "condition", value = "All Selected")
          updateTextInput(session, "background", value = "All Selected")
          
        }
        else if(input$View == "View By Selection Condition"){
          shinyjs::enable("cumulView")
          shinyjs::disable("classView")
          updateTextInput(session, "instructor", value = "All Selected")
          updateTextInput(session, "year", value = "All Selected")
          updateTextInput(session, "sample", value = "All Selected")
        }
      } else {
        shinyjs::disable("classView")
        shinyjs::disable("cumulView")
      }
    })
    
    #Handling behaviors for button selections
    observe({
      options <- c(as.character(mutation_data() %>% filter(condition == input$condition) %>% pull(background)))
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
        updateSelectizeInput(session, "year", choices = c("All Selected", unique(mutation_data()$year)), server=TRUE)
      } else {
        updateSelectizeInput(session, "year", choices = c("All Selected", as.character(mutation_data()[mutation_data()$instructor == input$instructor, "year"])), server=TRUE)
      }
    })
    
    observe({
      updateSelectInput(session, "instructor", choices = c('All Selected', unique(mutation_data()$instructor)))
    }) 
    
    observe({
      updateSelectInput(session, "sample", choices = c("All Selected", as.character(mutation_data() %>% filter(instructor==input$instructor) %>% filter(year==input$year) %>% pull(sample))))
    })
    
    observe({
      if(input$sample != "All Selected") {
        choices <- mutation_data()[mutation_data()$sample == input$sample, "GENE"]
      } else {
        choices <- filtered_data() %>% pull(GENE)
      }
      
      # Filter choices to include only those present in genes_info
      choices <- choices[choices %in% genes_info$GENE]
      
      choices <- sort(choices)
      
      updateSelectizeInput(session, "GENE", choices = as.character(choices), server=TRUE)
    })
    
    
    # Learn about Gene button within gene viewer
    sgdid <- reactiveValues(value = NULL)
    output$url <- renderUI({
      sgdid_values <- genes_info[genes_info$GENE == input$GENE,"SGDID"]
      sgdid$value <- sgdid_values
      url <- a("Learn about Gene",href=paste0(link, sgdid$value),class="btn btn-default", target='_blank')
      url
    })
    
    output$pdf_viewer <- renderUI({ tags$iframe(
      style="height:1000px;width:100%;scrolling=yes",
      src = "Black_box.pdf") }) 
    
    tabPanel("Ancestor Strain",
            uiOutput("pdf_viewer") )
    
    # to create loading message below: 
    loading_message <- "Loading..."
    # Calculate the number of empty spaces needed on each side
    total_spaces <- 160  # Total number of characters to occupy the line
    message_length <- nchar(loading_message)
    spaces_on_each_side <- floor((total_spaces - message_length) / 2)
    # Construct the string with spaces on each side of the loading message
    formatted_loading_message <- sprintf("%s%s%s%s",
                                          "\n\n\n\n\n\n\n\n",
                                          paste(rep(" ", spaces_on_each_side), collapse = ""),
                                          loading_message,
                                          paste(rep(" ", spaces_on_each_side), collapse = ""))
    
    # Define the desired order of categories
    desired_order <- c('chrM', 'chrXVI', 'chrXV', 'chrXIV', 'chrXIII', 'chrXII', 'chrXI', 'chrX', 'chrIX', 'chrVIII', 'chrVII', 'chrVI', 'chrV', 'chrIV', 'chrIII', 'chrII','chrI')
    # Convert category to a factor with the desired order
    chrom_info$CHROM <- factor(chrom_info$CHROM, levels = desired_order)

    # Define a mapping from chromosome names to numbers
    chromosome_mapping <- c(
      chrM = 1, chrXVI = 2, chrXV = 3, chrXIV = 4, chrXIII = 5, chrXII = 6,
      chrXI = 7, chrX = 8, chrIX = 9, chrVIII = 10, chrVII = 11, chrVI = 12,
      chrV = 13, chrIV = 14, chrIII = 15, chrII = 16, chrI = 17
    )
    
    # Create an empty dataframe to store the final results
    final_gene <- reactive({
      mutation_data_value <- filtered_data()
      # Merge the data frames based on the “REGION” column
      common_cols <- intersect(colnames(mutation_data_value), colnames(genes_info))
      mutation_data_value <- merge(mutation_data_value, genes_info, by = common_cols)
        
      # Iterate through unique genes
      final_gene_static <- mutation_data_value %>%
        group_by(GENE) %>%
        summarize(
          CHROM = first(CHROM),
          START = first(START),
          STOP = first(STOP),
          Counts = n(),
          chrom_as_num = first(chromosome_mapping[match(first(CHROM), names(chromosome_mapping))])
        ) %>%
        ungroup()

      return(final_gene_static)
    })
      
      
    output$chromPlot <- chromMapServer("chromPlot", final_gene, formatted_loading_message, chrom_info)
    output$varPieChart <- variantsServer("varPieChart", mutation_data, filtered_data)
    output$snpCountPlot <- snpCountServer("snpCountPlot", filtered_data)
    output$geneViewPlot <- geneViewServer("geneView", total_spaces, reactive(input$GENE), filtered_data, genes_info)
    
    
    observeEvent(input$append_btn, {
      new_csv_path <- input$new_csv$datapath
      
    }) 
  }

  shinyApp(ui, server)
}
