# If you are changing the 'masterfile' that the app is built off of (aka the initial database),
# This needs to be changed in two locations:
# mut_backend before the server is called, and mutation_data right after the server is called.


# List of required packages
required_packages <- c("devtools", "devtools", "shinythemes", "DBI", "RSQLite", 
                       "ggplot2","dplyr", "tidyr", "forcats", "ggrepel",
                       "purrr", "shinyjs", "viridis", "plotly")

# Check if packages are installed, and if not, install them
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

### %% Imports
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
library(stringr)

# loading in the VCF file to display initial choices, later turns into reactive val called mutation_data that includes manually updated data
# THIS SHOULD NOT BE CHANGED IN THE CODE. If the overall master shifts you can modify it here, 
# but mut_backend should not be written to anywere in the code it is designed to just fill in in the beginning, 
# right before the reactive frame gets created.
mut_backend <- read.csv("all_yEvo_vcf_spring2024.csv") 

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
               width = 3,
               fileInput("datafile", "Optional: Upload additional CSV File", accept = ".csv"),
               conditionalPanel(
    # Asks for instructor/year so we can modify user uploaded file so we can combine it
    # with our current data
                 condition = "output.filesUploaded",
                 textInput("inputted_instructor", "Who is your Instructor"),
                 textInput("inputted_year", "What is the current year"),
                 actionButton("submit_teach_year", "Submit Teacher and Year")
               ),
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
                 tabPanel("Chromosome Map", plotlyOutput("chromPlot",height = "600px"),verbatimTextOutput("info")),
                 tabPanel("Variant Pie Chart", plotlyOutput("varPieChart"), verbatimTextOutput("text")),
                 tabPanel("SNP Counts", plotlyOutput("snpCountPlot")),
                 tabPanel("Gene View", div("", style = "height: 10px;"), plotlyOutput("geneViewPlot", width = "600px"), verbatimTextOutput("gene"),
                         selectInput("GENE", "Gene", choices = c('')),
                         uiOutput("url")),
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
  #CHANGE MASTERFILE HERE IF NEEDED
  mutation_data <- reactiveVal(read.csv("all_yEvo_vcf_spring2024.csv")) 
  
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
    if(input$sample != "All Selected") {
      choices <- mutation_data()[mutation_data()$sample == input$sample, "GENE"]
    } else {
      choices <- filtered_data() %>% pull(GENE)
    }
    
    # Filter choices to include only those present in genes_info
    choices <- choices[choices %in% genes_info$GENE]
    
    choices <- sort(choices)
    
    updateSelectInput(session, "GENE", choices = as.character(choices))
  })
  
  
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
  
  # to create loading message below: 
  loading_message <- "Loading..."
  # Calculate the number of empty spaces needed on each side
  total_spaces <- 160  # Total number of characters to occupy the line
  message_length <- nchar(loading_message)
  spaces_on_each_side <- floor((total_spaces - message_length) / 2)
  # Construct the string with spaces on each side of the loading message
  formatted_loading_message <- sprintf("%s%s%s",
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
  
  output$chromPlot <- renderPlotly({
    validate(
      need(final_gene()$START, formatted_loading_message)
    )
    
    # Plotting
    final_gene_data <- final_gene()
    final_gene_data$CHROM <- factor(final_gene_data$CHROM, levels = levels(chrom_info$CHROM))
    
    scale_y_custom <- ggplot2::scale_y_continuous(
      breaks = rev(as.numeric(chrom_info$CHROM)),
      labels = rev(as.character(chrom_info$CHROM))
    )
    
    final_gene_singletons <- final_gene_data[final_gene_data$Counts == 1, ]
    if (any(final_gene_data$Counts >= 2)) {
      final_gene_multi_muts <- final_gene_data[final_gene_data$Counts >= 2, ]
    } 
    p <- ggplot() +
      geom_rect(data = chrom_info,
                aes(xmin = 0, xmax = length, ymin = as.numeric(factor(CHROM)) - .02, ymax = as.numeric(factor(CHROM)) + .02, text = CHROM
                    ),
                fill = 'lightblue', alpha = 1) +
      geom_rect(data = chrom_info,
                aes(xmin = 0, xmax = length, ymin = as.numeric(factor(CHROM)) - 0.2, ymax = as.numeric(factor(CHROM)) + 0.2,
                    text = CHROM),
                fill = 'lightblue', alpha = 1) +
      geom_rect(data = final_gene_singletons, aes(ymin = chrom_as_num - 0.4, # swapped ymin and ymax
                                                  ymax = chrom_as_num + 0.4,
                                                  xmin = START,
                                                  xmax = START + 8000,
                                                  text = paste0("Gene Name: ",GENE,"\n", "Independent Mutations: ",Counts)),
                fill = "white", alpha = 1, color = "black", size = 0.1) +
      scale_y_custom +
      scale_fill_gradient(low = "pink", high = "red4",) +
      labs(title = 'Location of mutations along chromosomes',
           y = 'Chromosome', # changed x-axis label to Chromosome
           x = 'Position along chromosome') # changed y-axis label to Length
    if (exists("final_gene_multi_muts")) {
      p <- p + geom_rect(data = final_gene_multi_muts, 
                         aes(ymin = chrom_as_num - 0.4, ymax = chrom_as_num + 0.4,
                             xmin = START, xmax = START + 8000,
                             text = paste0("Gene Name: ",GENE,"\n", "Independent Mutations: ",Counts), fill = Counts),
                         alpha = 1, color = "black", size = 0.1)
    }
      
    # Convert ggplot2 plot to plotly
    p <- ggplotly(p,tooltip = "text")
    # Add formatting
    p <- layout(
      p,
      plot_bgcolor = 'rgba(0,0,0,0)',   # Set plot background color to transparent
      paper_bgcolor = 'rgba(0,0,0,0)',  # Set paper background color to transparent
      xaxis = list(showgrid = FALSE),  # Remove x-axis gridlines
      yaxis = list(showgrid = FALSE)   # Remove y-axis gridlines
    )
  })
  
  # Render Pie Chart
  output$varPieChart <- renderPlotly({
    # Color vector for annotations 
    # just add the same number of colors as number of annotations
    # ex. if there are 10 unique annotations, put 10 unique colors
    # in this color vector
    color_vector <- c("#9edae5", "#17becf", "#dbdb8d", "#bcbd22", "#c7c7c7",
                      "#7f7f7f", "#f7b6d2", "#e377c2", "#c49c94", "#8c564b",
                      "#c5b0d5", "#9467bd", "#ff9896", "#d62728", "#98df8a",
                      "#2ca02c", "#ffbb78", "#ff7f0e", "#aec7e8", "#1f77b4")
    
    # Filter and get a vector of unique annotations (no duplicates)
    all_unique_anno <- mutation_data() %>%
      distinct(ANNOTATION) %>% pull(ANNOTATION)
    
    # sort the annotations
    all_unique_anno <- sort(all_unique_anno)
    
    # initialize a data frame which has 3 columns: Annotations, count of annotations, and percentage
    pie_df <- data.frame(
      ANNOTATION = all_unique_anno,
      count = numeric(length(all_unique_anno)),
      percent = numeric(length(all_unique_anno))
    )

    # filter the filtered data further:
    # sort the data alphabetically
    # count the number of annotations
    # get a percentage for each annotation
    pie_data <- filtered_data() %>%
      arrange(ANNOTATION) %>% 
      count(ANNOTATION) %>% 
      mutate(percent = n/sum(n) * 100)
    
    # Update pie_df with counts from filtered_data
    pie_df$count[match(pie_data$ANNOTATION, pie_df$ANNOTATION)] <- pie_data$n
    
    # Calculate percent
    pie_df$percent <- pie_df$count / sum(pie_df$count) * 100
    
    # Arrange the data alphabetically
    sorted_pie_df <- arrange(pie_df, ANNOTATION)
    
    p <- plot_ly(sorted_pie_df, labels = ~ANNOTATION, values = ~percent, type = 'pie', text = ~paste(ANNOTATION, ": ", round(percent, digits = 2), "%"),
                 hoverinfo = "text", outsidetextfont = list(size = 8), textinfo = "text", marker = list(colors = color_vector)) %>%
      layout(title = list(text = "Percentage of Variants by Type", x = 1, y = 0.95, xanchor = "right", yanchor = "top"),
             legend = list(x = 0.95, y = 0.1),
             plot_bgcolor = 'rgba(0,0,0,0)',
             paper_bgcolor = 'rgba(0,0,0,0)',
             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    
    p %>% layout(showlegend = TRUE, 
                 legend = list(font = list(size = 7)),
                 legend = list(font = list(size = 7)),
                 legend = list(font = list(size = 7)),
                 margin = list(l = 75, r = 75, b = 75, t = 75))
  })
  
  # Render SNP Plot
  output$snpCountPlot <- renderPlotly({
    # categorizing data into appropriate categories for plotting
    categorized_data <- filtered_data() %>%
      mutate(transition = paste(REF, " to ", ALT, sep = "")) %>%
      mutate(transition = if_else(nchar(transition) > 6,"Indel",transition)) %>%
      mutate(mutation_type = case_when(
        transition %in% c("A to G", "G to A", "C to T", "T to C") ~ "Transition",
        transition %in% c("A to T", "T to A", "C to G", "G to C", "A to C", "C to A", "T to G", "G to T") ~ "Transversion",
        TRUE ~ "Indel"
      ))
    
    num_categories <- length(unique(categorized_data$mutation_type)) # in case we ever change categories
    # Getting counts to display in hover (if not displaying custom tooltip then can use categorized data directly)
    summarized_data <- categorized_data %>%
      group_by(transition, mutation_type) %>%
      summarise(count = n(), .groups = 'drop')
    
    # plotting the data with coloring by categories
    p <- ggplot(summarized_data, aes(x = as.factor(transition), fill = mutation_type,)) +
      geom_bar(aes(y = count, text = paste(transition, '\nCount:', count,'\nMutation Type:', mutation_type)), stat = "identity") +
      theme_bw() +
      scale_fill_manual(values = viridis::viridis(num_categories, begin = 0.4, end = 1)) +
      labs(fill = "Mutation Type") +
      theme(
        legend.position = "right",
        strip.text.y.left = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=.5)
      ) +
      ggtitle("Single Nucleotide Changes") +
      xlab("SNP call")
    p <- ggplotly(p, tooltip = "text")
  })
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  # Showing please select gene message
  # Construct the string with spaces on each side of the loading message
  geneview_message <- "Please select a gene below"
  # Calculate the number of empty spaces needed on each side
  gene_message_length <- nchar(geneview_message)
  gene_spaces_on_each_side <- floor((total_spaces - gene_message_length) / 2)
  select_gene_message <- sprintf("%s%s%s",
                                       "\n\n\n\n\n\n\n\n",
                                       paste(rep(" ", gene_spaces_on_each_side), collapse = ""),
                                       geneview_message,
                                       paste(rep(" ", gene_spaces_on_each_side), collapse = ""))
  
  # Render gene view plot
  output$geneViewPlot <- renderPlotly({
    validate(
      need(input$GENE, select_gene_message)
    )
    xlength <- genes_info %>%
      filter(GENE==input$GENE) %>% pull(PROTEIN_LENGTH) %>% unique() %>% as.numeric()
    
    # TESTING
    gene_name <- input$GENE
    
    mutation_data_value <- filtered_data()
    
    # Merge the data frames based on the “REGION” column
    common_cols <- intersect(colnames(mutation_data_value), colnames(genes_info))
    mutation_data_value <- merge(mutation_data_value, genes_info, by = common_cols)
    
    mutation_data_value <- mutation_data_value[order(mutation_data_value$GENE),]
    
    # Using subset() function
    cur_gene <- subset(mutation_data_value, GENE == gene_name)
    
    # Iterate through unique genes
    count_proteins <- cur_gene %>%
      group_by(POS) %>%
      summarize(
        GENE = first(GENE),
        PROTEIN = first(PROTEIN),
        ANNOTATION = first(ANNOTATION),
        COUNTS = n(),
        Letter1 = substr(PROTEIN, 1, 1),  # Extract the first character
        Numbers = as.numeric(str_extract(PROTEIN, "[0-9]+")),  # Extract the numbers
        Letter2 = str_extract(PROTEIN, "[a-zA-Z]+$")
      ) %>%
      ungroup()
    
    count_proteins_same <- count_proteins %>%
      group_by(Numbers) %>%
      summarize(
        GENE = first(GENE),
        PROTEIN = paste((PROTEIN), collapse = ", "),
        ANNOTATION = first(ANNOTATION),
        Counts_diff_mutation = paste((COUNTS), collapse = ", "),
        Counts_tot = sum(COUNTS)
      ) %>%
      ungroup()
    count_proteins_same$PROTEIN <- sapply(count_proteins_same$PROTEIN, FUN = strsplit, split = ",", simplify = TRUE)
    count_proteins_same$Counts_diff_mutation <- sapply(count_proteins_same$Counts_diff_mutation, FUN = strsplit, split = ",", simplify = TRUE)
    
    combined_strings <- character(nrow(count_proteins_same))  # Pre-allocate character vector
    for (i in 1:nrow(count_proteins_same)) {
      current_protein <- (count_proteins_same$PROTEIN[i])
      current_counts <- count_proteins_same$Counts_diff_mutation[i]
      
      split_string <- strsplit(current_protein[[1]], ",")
      split_counts <- strsplit(current_counts[[1]], ",")
      
      if (length(split_string) > 1) {
        cur <- list()
        for (j in 1:length(split_string)) {
          current_protein = split_string[j]
          current_protein <- str_trim(current_protein)
          current_counts = split_counts[j]
          Letter1 <- substr(current_protein, 1, 1)  # Extract the first character
          Numbers <- as.numeric(str_extract(current_protein, "[0-9]+"))  # Extract the numbers
          Letter2 <- str_extract(current_protein, "[a-zA-Z]+$")
          print(current_protein)
          print(Letter1)
          cur <- c(cur, paste("Count ", Letter1, '->', Letter2, ": ", current_counts, "\n"))
        }
        combined_strings[i] <- sapply(cur, function(x) paste(x)) %>% paste(collapse = "")
      }
      else {
        Letter1 <- substr(current_protein, 1, 1)  # Extract the first character
        Numbers <- as.numeric(str_extract(current_protein, "[0-9]+"))  # Extract the numbers
        Letter2 <- str_extract(current_protein, "[a-zA-Z]+$")
        combined_strings[i] <- paste("Count ", Letter1, '->', Letter2, ": ", current_counts)
      }
    }
    count_proteins_same$combined <- combined_strings
    
    
#TESTING
    max_count <- max(count_proteins$COUNTS)
    
    p <- count_proteins_same %>%
      mutate(
        AA_WT = substr(PROTEIN, 1, 1),  # Extract the first character Amino Acid Wild Type
        AA_POS = as.numeric(str_extract(PROTEIN, "[0-9]+")),  # Extract Amino Acid Position
        AA_M = substr(PROTEIN, nchar(PROTEIN), nchar(PROTEIN)) # Amino Acid Mutation
      ) %>%
      mutate(ANNOTATION= gsub("'","",ANNOTATION)) %>%
      mutate(AA_POS = if_else(ANNOTATION=="5-upstream",-15,as.numeric(AA_POS))) %>%
      
      ggplot(aes(x = as.numeric(AA_POS), y = max_count + 4, 
                            text = ifelse(is.na(PROTEIN), 
                                   paste0('Non-coding Mutation\nAnnotation: ', ANNOTATION, '\nCount: ', Counts_diff_mutation), 
                                   paste0(combined, '\nPosition: ', AA_POS))))+
      geom_hline(yintercept=0, linetype=2,alpha=.2)+
      geom_segment(aes(x=0,xend=xlength,y=0,yend=0), size=15, color = "cornflowerblue") +
      geom_segment(aes(x=as.numeric(AA_POS),xend=as.numeric(AA_POS),y=0,yend=Counts_tot), color = "pink") +
      geom_point(aes(x=as.numeric(AA_POS), y=Counts_tot,color=ANNOTATION), size=2) +
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
            plot.title = element_text(hjust = 0.5),  # Adjust top margin for title
            axis.text.y = element_blank(),
            plot.margin = margin(20, 0, 0, 0) # Adjust top margin for space between title and plot
      ) + 
      xlab("Amino acid position") +
      theme(axis.text.x = element_text(size = 8)) +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) + 
      theme(plot.margin = margin(0, 0, 0, 0)) + 
      annotate(
        "text", x = 1, y = Inf, label = "Drag over mutations to see more",
        hjust = 0, vjust = 2, color = "black", size = 5
      ) + 
      guides(color = FALSE)
    
    ggplotly(p, tooltip="text")
  })
  
  
  observeEvent(input$append_btn, {
    new_csv_path <- input$new_csv$datapath
    
  }) 
}

shinyApp(ui, server)