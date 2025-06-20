initialize_dataframe <- function(input, mut_backend) {
  # Initialize a reactive variable for the dataframe
  # Function to read and append the uploaded data to the cumulative dataframe
  observeEvent(input$submit_teach_year, {
    file <- input$datafile
    data <- read.csv(file$datapath)
    required_columns <- c(
      "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "ANNOTATION", "REGION", "GENE", "PROTEIN", "seq_file", "background",
      "condition", "sample"
    )
    if (!is.null(file) && all(required_columns %in% colnames(data))) {
      # Add instructor and year columns
      data$instructor <- rep(input$inputted_instructor, nrow(data))
      data$year <- rep(input$inputted_year, nrow(data))
      # Rearrange columns
      data <- data[, c(
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "ANNOTATION", "REGION", "GENE", "PROTEIN", "seq_file",
        "background", "condition", "instructor", "year", "sample"
      )]
      df <- rbind(mut_backend, data)
      mutation_data(df)
    } else {
      # Handle the case where required columns are missing
      print("Some required columns are missing.")
      mutation_data(mut_backend)
    }
  })
}

set_visualization_settings <- function(input, session, genes_info, mutation_data, filtered_data) {
  # Display settings
  observe({
    if (!is.null(input$View)) {
      if (input$View == "View By Class") {
        shinyjs::disable("cumulView")
        shinyjs::enable("classView")
        updateTextInput(session, "condition", value = "All Selected")
        updateTextInput(session, "background", value = "All Selected")
      } else if (input$View == "View By Selection Condition") {
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
  
  # Handling behaviors for button selections
  observe({
    options <- c(as.character(mutation_data() %>%
                                filter(condition == input$condition) %>% pull(background)))
    if (length(unique(options)) != 1) {
      updateSelectInput(session, "background",
                        choices = c("All Selected", options)
      )
      shinyjs::enable("background")
    } else {
      updateSelectInput(session, "background", choices = c(options))
      shinyjs::disable("background")
    }
  })
  # making the background not a button before selecting a condition
  observe({
    if (input$condition == "All Selected") {
      shinyjs::disable("background")
    }
  })
  
  observe({
    updateSelectInput(session, "instructor", choices = c(
      "All Selected",
      unique(mutation_data()$instructor)
    ))
  })
  
  observe({
    # Make sure the data exists before we do anything
    req(mutation_data())
    
    # Figure out which years to show
    if (input$instructor == "All Selected") {
      years <- sort(unique(as.character(mutation_data()$year)))
    } else {
      years <- sort(unique(as.character(
        mutation_data()$year[
          mutation_data()$instructor == input$instructor
        ]
      )))
    }
    
    # Update the dropdown menu with clean year choices
    updateSelectizeInput(
      session,
      "year",
      choices = c("All Selected", years),
      selected = "All Selected",
      server = TRUE
    )
  })
  
  
  observe({
    updateSelectInput(session, "instructor", choices = c(
      "All Selected",
      unique(mutation_data()$instructor)
    ))
  })
  
  observe({
    updateSelectInput(session, "sample", choices = c(
      "All Selected",
      as.character(mutation_data() %>%
                     filter(instructor == input$instructor) %>%
                     filter(year == input$year) %>%
                     pull(sample))
    ))
  })
  
  #gene view draopdown menu
  observe({
    if (input$sample != "All Selected") {
      choices <- mutation_data()[
        mutation_data()$sample == input$sample,
        "GENE"
      ]
    } else {
      choices <- filtered_data() %>% pull(GENE)
    }
    
    # Filter choices to include only those present in genes_info
    choices <- choices[choices %in% genes_info$GENE]
    
    choices <- sort(choices)
    
    updateSelectizeInput(session, "GENE",
                         choices = as.character(choices),
                         server = TRUE,
                         options = list(maxOptions = length(choices))
    )
  })
}
