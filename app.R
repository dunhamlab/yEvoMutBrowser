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
library(scales)

# loading in the VCF file to display initial choices, later turns into reactive val called mutation_data that includes manually updated data
# THIS SHOULD NOT BE CHANGED IN THE CODE. If the overall master shifts you can modify it here, 
# but mut_backend should not be written to anywere in the code it is designed to just fill in in the beginning, 
# right before the reactive frame gets created.
mut_backend <- read.csv("all_yEvo_vcf_spring2025.csv") 

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
                         selectInput("GENE", "Gene", choices = NULL),
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
  mutation_data <- reactiveVal(read.csv("all_yEvo_vcf_spring2025.csv")) 
  
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
  
  tabPanel("Background",
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
    
    
  output$chromPlot <- renderPlotly({
    validate(
      need(final_gene()$START, formatted_loading_message)
    )
    
    # Preparing data
    final_gene_data <- final_gene()
    final_gene_data$CHROM <- factor(final_gene_data$CHROM, levels = levels(chrom_info$CHROM))
    
    # Splitting data so that the colors don't get messed up with the single counts
    final_gene_singletons <- final_gene_data %>% filter(Counts == 1)
    final_gene_multi_muts <- final_gene_data %>% filter(Counts >= 2)
    
    # Add numerical position for graphing chromosomes
    final_gene_singletons$chrom_num <- as.numeric(final_gene_singletons$CHROM) - 1
    final_gene_multi_muts$chrom_num <- as.numeric(final_gene_multi_muts$CHROM) - 1
    
    # Creating chromosome bars
    fig <- plot_ly(
      data = chrom_info,
      x = ~length, 
      y = ~CHROM,
      type = "bar",
      orientation = 'h',
      marker = list(color = 'lightblue'),
      name = 'Chromosomes',
      text = ~CHROM,
      textposition = 'none',
      hoverinfo = 'text',
      showlegend = FALSE
    )
    fig <- fig %>%
      layout(
        bargap = 0.3  # giving a little more space between chrom bars
      )
    
    
    # Shapes list for mutation markers
    shapes_list <- list()
    
    # Adding singleton mutation rectangles
    if (nrow(final_gene_singletons) > 0) {
      for (i in 1:nrow(final_gene_singletons)) {
        shapes_list[[length(shapes_list) + 1]] <- list(
          type = "rect",
          x0 = final_gene_singletons$START[i], 
          x1 = final_gene_singletons$START[i] + 8000,
          y0 = final_gene_singletons$chrom_num[i] - 0.4,
          y1 = final_gene_singletons$chrom_num[i] + 0.4,
          fillcolor = "white",
          line = list(color = "black", width = 0.1),
          layer = "above"
        )
      }
      
      # Adding invisible points for singleton hover functionality
      fig <- fig %>% 
        add_trace(
          data = final_gene_singletons,
          x = ~START + 4000,
          y = ~CHROM,
          type = "scatter",
          mode = "markers",
          marker = list(size = 10, opacity = 0,color = "white"),
          name = 'Single Mutations',
          text = ~paste0("Gene Name: ", GENE, "<br>Independent Mutations: ", Counts),
          hoverinfo = 'text'
        )
    }
    
    # Add multi-mutation rectangles
    if (nrow(final_gene_multi_muts) > 0) {
      # Create color scale
      max_count <- max(final_gene_multi_muts$Counts)
      # Create a vector to store colors for each multi mutation
      multi_mut_colors <- character(nrow(final_gene_multi_muts))
      
      for (i in 1:nrow(final_gene_multi_muts)) {
        # Calculate color based on count
        count_ratio <- (final_gene_multi_muts$Counts[i] - 2) / (max_count - 2)
        if (is.na(count_ratio) || count_ratio < 0) count_ratio <- 0
        
        # Linear interpolation between pink and dark red to match default color interpolation on points
        r <- 255 - count_ratio * (255 - 139)  # pink(255) to red4(139)
        g <- 192 - count_ratio * 192          # pink(192) to red4(0)
        b <- 203 - count_ratio * 203          # pink(203) to red4(0)
        
        hex_color <- rgb(r/255, g/255, b/255)
        multi_mut_colors[i] <- hex_color
        
        shapes_list[[length(shapes_list) + 1]] <- list(
          type = "rect",
          x0 = final_gene_multi_muts$START[i], 
          x1 = final_gene_multi_muts$START[i] + 8000,
          y0 = final_gene_multi_muts$chrom_num[i] - 0.4,
          y1 = final_gene_multi_muts$chrom_num[i] + 0.4,
          fillcolor = hex_color,
          line = list(color = "black", width = 0.1),
          layer = "above"
        )
      }
      final_gene_multi_muts$color <- multi_mut_colors
      
      # Adding invisible points for multi mutation hover functionality
      fig <- fig %>% 
        add_trace(
          data = final_gene_multi_muts,
          x = ~START + 4000,
          y = ~CHROM,
          type = "scatter",
          mode = "markers",
          marker = list(
            color = ~Counts,  # assuming you want to map 'Counts' to color
            colorscale = list(
              list(0, "#FFC0CB"),   
              list(1, "#8B0000")    
            ),
            size = 10,
            opacity = 0,
            colorbar = list(title = "Mutation Count")
            ),
          name = 'Multiple Mutations',
          text = ~paste0("Gene Name: ", GENE, "<br>Independent Mutations: ", Counts),
          hoverinfo = 'text'
        )
    }
  
    
    # Addding shapes to graph
    fig <- fig %>% layout(
      title = list(
        text = 'Location of mutations along chromosomes',
        x = 0.5,
        xanchor = 'center',
        y = 0.95,
        yanchor = 'top',
        font = list(size = 20)
      ),
      xaxis = list(
        title = 'Position along chromosome',
        showgrid = FALSE
      ),
      yaxis = list(
        title = 'Chromosome',
        showgrid = FALSE,
        autorange = "reversed"
      ),
      shapes = shapes_list,
      plot_bgcolor = 'rgba(0,0,0,0)',
      paper_bgcolor = 'rgba(0,0,0,0)'
    )
    
    fig
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
  
  output$snpCountPlot <- renderPlotly({
    # Categorizing data into appropriate categories for plotting
    categorized_data <- filtered_data() %>%
      mutate(transition = paste(REF, " to ", ALT, sep = "")) %>%
      mutate(transition = if_else(nchar(transition) > 6, "Indel", transition)) %>%
      mutate(mutation_type = case_when(
        transition %in% c("A to G", "G to A", "C to T", "T to C") ~ "Transition",
        transition %in% c("A to T", "T to A", "C to G", "G to C", "A to C", "C to A", "T to G", "G to T") ~ "Transversion",
        TRUE ~ "Indel")) %>%
      mutate(combined_group = case_when(
        transition %in% c("A to G", "T to C") ~ "A to G",
        transition %in% c("C to T", "G to A") ~ "C to T",
        transition %in% c("A to T", "T to A") ~ "A to T",
        transition %in% c("C to G", "G to C") ~ "C to G",
        transition %in% c("A to C", "T to G") ~ "A to C",
        transition %in% c("G to T", "C to A") ~ "C to A",
        TRUE ~ "Indel"
      ))
    
    # Summarize the data by the new group variable
    summarized_data <- categorized_data %>%
      group_by(combined_group, mutation_type) %>%
      summarise(count = n(), .groups = 'drop')
    
    # Custom colors for mutation types
    custom_colors <- c(
      "Transition" = "#0072B2", # Blue
      "Transversion" = "#CC79A7", # Pink
      "Indel" = "#009E73" # Green        
    )
    desired_order <- c('Indel', 'A to G', 'C to T', 'A to T', 'C to G', 'A to C', 'C to A')
    
    summarized_data$combined_group <- factor(summarized_data$combined_group, levels = desired_order)
    
    # Plotting the data with coloring by categories
    p <- ggplot(summarized_data, aes(x = (combined_group), y = count, fill = mutation_type)) +
      aes(text = paste(combined_group, '\nCount:', count, '\nMutation Type:', mutation_type)) +
      geom_bar(position = "dodge", stat = "identity") +
      theme_bw() +
      scale_fill_manual(values = custom_colors) +
      labs(fill = "Mutation Type") +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) # Adjust angle for better visibility
      ) +
      ggtitle("Single Nucleotide Changes") +
      xlab("Mutation Type and SNP Call")
    
    p <- ggplotly(p, tooltip = "text")
  })
  
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  # Showing please select gene message
  # Construct the string with spaces on each side of the loading message
  geneview_message <- "Please select a gene below"
  # Calculate the number of empty spaces needed on each side
  gene_message_length <- nchar(geneview_message)
  gene_spaces_on_each_side <- floor((total_spaces - gene_message_length) / 2)
  select_gene_message <- sprintf("%s%s%s%s",
                                       "\n\n\n\n\n\n\n\n",
                                       paste(rep(" ", gene_spaces_on_each_side), collapse = ""),
                                       geneview_message,
                                       paste(rep(" ", gene_spaces_on_each_side), collapse = ""))
  
  # Render gene view plot
  output$geneViewPlot <- renderPlotly({
    validate(
      need(input$GENE, select_gene_message)
    )
    
    mutation_data_value <- filtered_data()
    
    # Merge the data frames based on the common columns
    common_cols <- intersect(colnames(mutation_data_value), colnames(genes_info))
    mutation_data_value <- merge(mutation_data_value, genes_info, by = common_cols)
    mutation_data_value <- mutation_data_value[order(mutation_data_value$GENE),]
    
    # Filter data for the specific gene
    cur_gene <- filter(mutation_data_value, GENE == input$GENE)
    
    # Pattern for extracting the second part of protein
    pattern <- "(?<=\\d)([A-Za-z]|\\*|indel)$|([A-Za-z]|\\*)$"
    
    all_annotations <- c("missense", "nonsense", "5'-upstream", "indel-frameshift","indel-inframe","synonymous")
    # IF ADDING NEW ANNOTATIONS - DON'T FORGET TO ADD BOTH HERE AND IN annotation_colors BELOW
    
    # Group and summarize protein counts
    count_proteins <- cur_gene %>%
      mutate(indel = nchar(ALT) - nchar(REF)) %>%  # Calculate indel difference
      group_by(POS) %>%
      summarize(
        GENE = first(GENE),
        PROTEIN = first(PROTEIN),
        ANNOTATION = first(ANNOTATION),
        COUNTS = n(),
        Letter1 = substr(PROTEIN, 1, 1),  # Extract the first character
        Numbers = as.numeric(str_extract(PROTEIN, "[0-9]+")),  # Extract the numbers
        Letter2 = str_extract(PROTEIN, pattern),
        indel = first(indel)
      ) %>%
      ungroup()
    
    # Group by position numbers and summarize
    count_proteins_same <- count_proteins %>%
      group_by(Numbers) %>%
      summarize(
        GENE = first(GENE),
        PROTEIN = list(PROTEIN),
        ANNOTATION = first(ANNOTATION),
        Counts_diff_mutation = list(COUNTS),
        Counts_tot = sum(COUNTS),
        indel = first(indel)
      ) %>%
      ungroup()
    
    # Combine protein and count strings
    count_proteins_same <- count_proteins_same %>%
      mutate(
        combined = map2_chr(PROTEIN, Counts_diff_mutation, function(prot, counts) {
          prot_list <- str_split(prot, ", ") %>% unlist()
          counts_list <- str_split(counts, ", ") %>% unlist()
          
          combined_strings <- map2_chr(prot_list, counts_list, function(p, c) {
            Letter1 <- substr(p, 1, 1)
            Numbers <- str_extract(p, "[0-9]+") %>% as.numeric()
            Letter2 <- str_extract(p, pattern)
            paste("Count", Letter1, '->', Letter2, ":", c, "\n", sep = " ")
          })
          paste(combined_strings, collapse = "")
        }) 
        
      )%>%
      mutate(
        PROTEIN = as.character(PROTEIN),
        AA_WT = substr(PROTEIN, 1, 1),  # Extract the first character Amino Acid Wild Type
        AA_POS = if_else(ANNOTATION == "5'-upstream", -15, as.numeric(str_extract(PROTEIN, "[0-9]+"))),  # Extract Amino Acid Position
        AA_M = substr(PROTEIN, nchar(PROTEIN), nchar(PROTEIN)),  #Amino Acid Mutation
        ANNOTATION = factor(ANNOTATION, levels = all_annotations) 
      )
    
    # Fixed colors for the different annotations - IF ADDING NEW ANNOTATIONS, DON'T FORGET TO ADD TO all_annotations ABOVE
    annotation_colors <- c(
      "missense" = "#800080",
      "nonsense" = "#61D04F",
      "5'-upstream" = "#F5C710",
      "indel-frameshift" = "#DF536B",
      "indel-inframe" = "#FFB6C1",
      "synonymous" = "cornflowerblue"
    )
    
    # Generating ranges for the plot size
    xmax <- genes_info %>%
      filter(GENE==input$GENE) %>% pull(PROTEIN_LENGTH) %>% unique() %>% as.numeric()
    ranges$x <- c(-50, xmax + 50)
    
    ranges$y <- c(0, max(count_proteins_same$Counts_tot) + ifelse(max(count_proteins_same$Counts_tot) < 5, 4, 1))

    p <- count_proteins_same %>%
      ggplot(
        aes(
        x = AA_POS,
        y = Counts_tot,
      )
      ) +
      # dummy geom_point to get legend showing all_annotations regardless of what data is displayed
      geom_point(data = data.frame(
          PROTEIN = rep(NA_character_, length(all_annotations)),
          AA_WT = rep(NA_character_, length(all_annotations)),
          AA_POS = rep(NA_real_, length(all_annotations)),
          AA_M = rep(NA_character_, length(all_annotations)),
          ANNOTATION = factor(all_annotations, levels = all_annotations),
          Counts_diff_mutation = I(rep(list(numeric(0)), length(all_annotations))),
          Counts_tot = rep(NA_integer_, length(all_annotations)),
          combined = rep(NA_character_, length(all_annotations))
      ), aes(y = 0, color = ANNOTATION, text = NULL), size = 2) +
      
      geom_hline(yintercept = 0, linetype = 2, alpha = .2,aes(text = NULL)) +
      geom_segment(aes(x = 0, xend = xmax, y = 0, yend = 0, text = NULL), linewidth = 15, color = "cornflowerblue") +
      geom_segment(aes(x = AA_POS, xend = AA_POS, y = 0, yend = Counts_tot, text = NULL), color = "gray") +
      geom_point( size = 2) +
      aes(x = AA_POS, 
           y = Counts_tot, 
           color = ANNOTATION,
           text = ifelse(
             grepl("indel", ANNOTATION),  # Check if ANNOTATION contains "indel"
             paste0("Indel\n", abs(indel), " base ",ifelse(indel > 0, "insertion", "deletion"), "\nCount: ", Counts_diff_mutation, "\nPosition: ", AA_POS),  # Text for indel annotations
             ifelse(
               is.na(PROTEIN), 
               paste0(ANNOTATION, '\nCount: ', Counts_diff_mutation, '\nPosition: ', -abs(unique(cur_gene$POS) - unique(cur_gene$START))),
               paste0(combined, '\nPosition: ', AA_POS)
             )
           )
      )+
      ggtitle(as.character(input$GENE)) +
      theme_classic() +
      theme(
        axis.text.x = element_text(),
        axis.title.y = element_text(),
        axis.text.y = element_text(),
        axis.ticks.y = element_line(),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(20, 20, 0, 20)
      ) +
      xlab("Amino acid position") +
      ylab("Mutation Count") +
      scale_y_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
      scale_color_manual(values = annotation_colors) +  # Manual color scale
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +      annotate("text", x = 1, y = Inf, label = "Drag over mutations to see more", hjust = 0, vjust = 2, color = "black", size = 5) +
      guides(color = guide_legend(title = "Annotation"))
    
    ggplotly(p, tooltip = "text")
  })
  
  
  observeEvent(input$append_btn, {
    new_csv_path <- input$new_csv$datapath
    
  }) 
}

shinyApp(ui, server)
