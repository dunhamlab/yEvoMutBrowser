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
      selectInput(NS(id, "instructor"), "Instructor", choices = ""),
      selectInput(NS(id, "year"), "Year", choices = ""),
      # PLEASE NOTE: "sample" is called "Sample Name" within ui to make it
      # easier to understand, but the official name in the df is sample
      selectInput(NS(id, "sample"), "Sample Name", choices = ""),
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

selection_panel_server <- function(id, filtered_data, mutation_data, mut_backend) {
  moduleServer(id, function(input, output, session) {
    # Initialize a reactive variable for the dataframe
    # Function to read and append the uploaded data to the cumulative dataframe
    observeEvent(input$submit_teach_year, {
      file <- input$datafile
      data <- read.csv(file$datapath)
      required_columns <- c(
        "CHROM", "POS", "REF", "ALT", "ANNOTATION", "REGION",
        "GENE", "PROTEIN", "background", "condition", "sample"
      )
      if (!is.null(file) && all(required_columns %in% colnames(data))) {
        # Add instructor and year columns
        data$instructor <- rep(input$inputted_instructor, nrow(data))
        data$year <- rep(input$inputted_year, nrow(data))
        data$ID <- rep("NA", nrow(data))
        data$QUAL <- rep("NA", nrow(data))
        data$FILTER <- rep("NA", nrow(data))
        data$INFO <- rep("NA", nrow(data))
        data$seq_file <- rep("NA", nrow(data))
        data <- subset(data, select=c(
          "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
          "ANNOTATION", "REGION", "GENE", "PROTEIN", "seq_file",
          "background", "condition", "instructor", "year", "sample"
        ))
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


    # download button functionality
    output$downloadBtn <- downloadHandler(
      filename = function() {
        # Set the filename for the downloaded file
        if (is.null(input$View) ||
          (input$instructor == "All Selected" &&
            input$condition == "All Selected")) {
          "master_table.csv"
        } else if (input$View == "View By Class") {
          selected_instructor <- input$instructor
          selected_year <- input$year
          selected_sample <- input$sample
          if (selected_sample != "All Selected") {
            paste0(
              selected_instructor, "_", selected_year, "_",
              selected_sample, ".csv"
            )
          } else if (selected_year != "All Selected") {
            paste0(selected_instructor, "_", selected_year, ".csv")
          } else if (selected_instructor != "All Selected") {
            paste0(selected_instructor, ".csv")
          }
        } else if (input$View == "View By Selection Condition") {
          if (input$background != "All Selected") {
            paste0(input$condition, "_", input$background, ".csv")
          } else {
            paste0(input$condition, ".csv")
          }
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
      if (input$instructor != "All Selected") {
        data <- data %>% filter(instructor == input$instructor)
      }
      if (input$year != "All Selected") {
        data <- data %>% filter(year == input$year)
      }
      if (input$sample != "All Selected") {
        data <- data %>% filter(sample == input$sample)
      }
      if (input$condition != "All Selected") {
        data <- data %>% filter(condition == input$condition)
      }
      if (input$background != "All Selected") {
        data <- data %>% filter(background == input$background)
      }
      data
    })

    server_outputs <- list(
      selected_instructor = reactive(input$instructor),
      selected_year = reactive(input$year),
      selected_sample = reactive(input$sample),
      selected_condition = reactive(input$condition),
      selected_background = reactive(input$background),
      filtered_data = new_filtered_data
    )

    return(server_outputs)
  })
}