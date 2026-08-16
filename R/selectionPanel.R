selection_panel_ui <- function(id, mut_backend) {
  sidebarPanel(
    width = 3,
    fileInput(NS(id, "datafile"), "Optional: Upload additional CSV File",
              accept = ".csv"
    ),
    conditionalPanel(
      # Asks for instructor/year so we can modify user uploaded file so we can
      # combine it with our current data
      condition = "output.filesUploaded",
      textInput(NS(id, "inputted_instructor"), "Who is your Instructor"),
      textInput(NS(id, "inputted_year"), "What is the current year"),
      actionButton(NS(id, "submit_teach_year"), "Submit Teacher and Year"),
      ns = NS(id)
    ),
    radioButtons(NS(id, "View"), "Select an option:",
                 choices = c("View By Class", "View By Selection Condition"),
                 selected = character(0)
    ),
    div("", style = "height: 10px;"), # Create a 10px vertical space
    
    # create download button here
    downloadButton(NS(id, "downloadBtn"), "Download the following data"),
    div("", style = "height: 10px;"), # Create a 10px vertical space
    
    # Only shows on condition observeEvent
    conditionalPanel(
      # links condition to button via button key
      condition = "input.uploadData || input.classView ||
      output.selectedClassView",
      # NOTE: instructor/year/sample are multi-select. "All Selected" is NOT
      # a real choice here -- it's shown as placeholder text (grayed out)
      # whenever nothing is picked, which is treated as "no filter" in the
      # server logic. This makes "pick specific values" and "all selected"
      # mutually exclusive by construction: there's no chip for "All
      # Selected" that could be picked alongside real values.
      selectizeInput(NS(id, "instructor"), "Instructor", choices = "",
                     multiple = TRUE,
                     options = list(placeholder = "All Selected")),
      selectizeInput(NS(id, "year"), "Year", choices = "",
                     multiple = TRUE,
                     options = list(placeholder = "All Selected")),
      # PLEASE NOTE: "sample" is called "Sample Name" within ui to make it
      # easier to understand, but the official name in the df is sample
      selectizeInput(NS(id, "sample"), "Sample Name", choices = "",
                     multiple = TRUE,
                     options = list(placeholder = "All Selected")),
      ns = NS(id)
    ),
    conditionalPanel(
      condition = "input.cumulView || output.selectedCumulView",
      selectInput(NS(id, "condition"), "Condition", choices = c(
        "All Selected",
        mut_backend %>%
          count(condition) %>%
          pull(condition)
      )),
      selectInput(NS(id, "background"), "Ancestor Strain", choices = ""),
      ns = NS(id)),
  )
}

region2gene_name <- function(gene_region, gene_info) {
  gene_info[gene_info$REGION == gene_region, "GENE"][1]
}

# Treats a multi-select input as "unfiltered" if nothing is picked yet.
# (For instructor/year/sample, "All Selected" is placeholder text rather
# than a real choice, so "nothing picked" IS the all-selected state.)
is_all_selected <- function(x) {
  is.null(x) || length(x) == 0
}

selection_panel_server <- function(id, filtered_data, mutation_data, mut_backend, gene_info) {
  moduleServer(id, function(input, output, session) {
    # Initialize a reactive variable for the dataframe
    # Function to read and append the uploaded data to the cumulative dataframe
    observeEvent(input$submit_teach_year, {
      file <- input$datafile
      data <- read.csv(file$datapath)
      required_columns <- c(
        "CHROM", "POS", "REF", "ALT", "ANNOTATION", "REGION", "PROTEIN", "background", "condition", "sample"
      )
      if (!is.null(file) && all(required_columns %in% colnames(data))) {
        # Add instructor and year columns
        data$instructor <- rep(input$inputted_instructor, nrow(data))
        data$year <- rep(input$inputted_year, nrow(data))
        data$INFO <- rep("NA", nrow(data))
        data$seq_file <- rep("NA", nrow(data))
        
        if (!('GENE' %in% data)) {
          data$GENE <- sapply(data$REGION, function(region_name) {region2gene_name(region_name, gene_info)})
        }
        data <- subset(data, select=c(
          "CHROM", "POS", "REF", "ALT", "INFO",
          "ANNOTATION", "REGION", "GENE", "PROTEIN", "seq_file",
          "background", "condition", "instructor", "year", "sample"
        ))
        # Rearrange columns
        data <- data[, c(
          "CHROM", "POS", "REF", "ALT", "INFO",
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
    
    # Display settings
    observe({
      if (!is.null(input$View)) {
        if (input$View == "View By Class") {
          shinyjs::disable("cumulView")
          shinyjs::enable("classView")
          updateSelectInput(session, "condition", selected = "All Selected")
          updateSelectInput(session, "background", selected = "All Selected")
        } else if (input$View == "View By Selection Condition") {
          shinyjs::enable("cumulView")
          shinyjs::disable("classView")
          # Reset the class-view multi-selects back to empty, which shows
          # the "All Selected" placeholder again.
          updateSelectizeInput(session, "instructor", selected = character(0))
          updateSelectizeInput(session, "year", selected = character(0))
          updateSelectizeInput(session, "sample", selected = character(0))
        }
      } else {
        shinyjs::disable("classView")
        shinyjs::disable("cumulView")
      }
    })
    
    # Handling behaviors for button selections
    observe({
      options <- as.character(mutation_data() %>%
                                pull(background))
      if (input$condition != "All Selected") {
        options <- as.character(mutation_data() %>%
                                  filter(condition == input$condition) %>%
                                  pull(background))
      }
      if (length(unique(options)) != 1) {
        updateSelectInput(session, "background",
                          choices = c("All Selected", options)
        )
        shinyjs::enable("background")
      } else {
        updateSelectInput(session, "background", choices = options)
        shinyjs::disable("background")
      }
    })
    
    observe({
      updateSelectizeInput(session, "instructor", choices = unique(mutation_data()$instructor))
    })
    
    observe({
      # Make sure the data exists before we do anything
      req(mutation_data())
      
      # Figure out which years to show, based on any/all instructors selected
      if (is_all_selected(input$instructor)) {
        years <- sort(unique(as.character(mutation_data()$year)))
      } else {
        years <- sort(unique(as.character(
          mutation_data()$year[
            mutation_data()$instructor %in% input$instructor
          ]
        )))
      }
      
      # Update the dropdown menu with clean year choices. No "selected" is
      # forced here, so the "All Selected" placeholder still shows until the
      # user actually picks something.
      updateSelectizeInput(
        session,
        "year",
        choices = years,
        server = TRUE
      )
    })
    
    observe({
      # Sample choices narrow based on whichever instructors/years are picked
      # (any number of them, or nothing/placeholder for no filter).
      updateSelectizeInput(session, "sample", choices = as.character(
        mutation_data() %>%
          filter(is_all_selected(input$instructor) | instructor %in% input$instructor) %>%
          filter(is_all_selected(input$year) | year %in% input$year) %>%
          pull(sample)
      ))
    })
    
    
    # download button functionality
    output$downloadBtn <- downloadHandler(
      filename = function() {
        # Build a compact label for a multi-select field: "All" if unfiltered,
        # otherwise the picked values joined with "-".
        fmt <- function(x) {
          if (is_all_selected(x)) "All" else paste(x, collapse = "-")
        }
        
        if (is.null(input$View)) {
          "master_table.csv"
        } else if (input$View == "View By Class") {
          paste0(
            fmt(input$instructor), "_",
            fmt(input$year), "_",
            fmt(input$sample), ".csv"
          )
        } else if (input$View == "View By Selection Condition") {
          if (input$background != "All Selected") {
            paste0(input$condition, "_", input$background, ".csv")
          } else {
            paste0(input$condition, ".csv")
          }
        } else {
          "master_table.csv"
        }
      },
      content = function(file) {
        # Write the data to a CSV file
        write.csv(filtered_data(), file)
      }
    )
    
    # storing a bool to see if a file has been uploaded
    # if a file has be uploaded, using the condition that if
    # output$filesUploaded
    # is true, we can auto-open the class view
    output$filesUploaded <- reactive({
      val <- !(is.null(input$datafile))
    })
    outputOptions(output, "filesUploaded", suspendWhenHidden = FALSE)
    
    output$selectedClassView <- reactive({
      if (!is.null(input$View)) {
        value <- (input$View == "View By Class")
      }
    })
    outputOptions(output, "selectedClassView", suspendWhenHidden = FALSE)
    
    output$selectedCumulView <- reactive({
      if (!is.null(input$View)) {
        value <- (input$View == "View By Selection Condition")
      }
    })
    outputOptions(output, "selectedCumulView", suspendWhenHidden = FALSE)
    
    # filtering dataframe based on menu selection, most things from here on out
    # should be based on filtered_data()
    new_filtered_data <- reactive({
      data <- mutation_data()
      # filtering based on selections if NOT all selected
      # (instructor/year/sample can each hold multiple picked values now;
      # an empty selection means "All Selected" / no filter)
      if (!is_all_selected(input$instructor)) {
        data <- data %>% filter(instructor %in% input$instructor)
      }
      if (!is_all_selected(input$year)) {
        data <- data %>% filter(year %in% input$year)
      }
      if (!is_all_selected(input$sample)) {
        data <- data %>% filter(sample %in% input$sample)
      }
      if (input$condition != "All Selected") {
        data <- data %>% filter(condition == input$condition)
      }
      if (input$background != "All Selected") {
        data <- data %>% filter(background == input$background)
      }
      data
    })
    
    # TRUE once the user has picked a View. For View By Class, an empty
    # multi-select is a valid, deliberate "All Selected" state (shown via
    # placeholder text), so it counts as complete just like a real pick.
    form_complete <- reactive({
      chosen_single <- function(x) !is.null(x) && nzchar(x) && x != "All Selected"
      
      v <- input$View
      if (is.null(v) || !nzchar(v)) return(FALSE)
      if (v == "View By Class") {
        TRUE
      } else if (v == "View By Selection Condition") {
        chosen_single(input$condition) && chosen_single(input$background)
      } else {
        FALSE
      }
    })
    
    server_outputs <- list(
      selected_instructor = reactive(input$instructor),
      selected_year = reactive(input$year),
      selected_sample = reactive(input$sample),
      selected_condition = reactive(input$condition),
      selected_background = reactive(input$background),
      filtered_data = new_filtered_data,
      form_complete = form_complete
    )
    
    return(server_outputs)
  })
}