classroom_dropdown_text <- "View yEvo Classroom data"
fhcc_sep_yevo_module_dropdown_text <- "View FHCC SEP yEvo Module data"

selection_panel_ui <- function(id, mut_backend) {
  all_year_options <- mut_backend$year %>% unique
  all_instructor_options <- mut_backend$instructor %>% unique
  all_condition_options <- mut_backend$condition %>% unique
  all_background_options <- mut_backend$background %>% unique
  all_sample_options <- mut_backend$sample %>% unique

  sidebarPanel(
    width = 3,
    # fileInput(NS(id, "datafile"), "Optional: Upload additional CSV File",
    #           accept = ".csv"
    # ),
    h4("Optional: Upload additional VCF data"),
    actionButton(NS(id, "uploadData"), "Upload VCF data"),
    # conditionalPanel(
    #   # Asks for instructor/year so we can modify user uploaded file so we can
    #   # combine it with our current data
    #   condition = "output.filesUploaded",
    #   textInput(NS(id, "inputted_instructor"), "Who is your Instructor"),
    #   textInput(NS(id, "inputted_year"), "What is the current year"),
    #   actionButton(NS(id, "submit_teach_year"), "Submit Teacher and Year"),
    #   ns = NS(id)
    # ),
    h4("Select data to view"),
    radioButtons(NS(id, "View"), "Select an option:",
                 choices = c(classroom_dropdown_text), # fhcc_sep_yevo_module_dropdown_text),
                 selected = classroom_dropdown_text
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
      pickerInput(
         NS(id, "instructor"),
         "Instructor",
         choices = all_instructor_options,
         selected = all_instructor_options,
         multiple = TRUE,
         options = pickerOptions(
           actionsBox = TRUE,
           selectedTextFormat = "count > 3"
         )
      ),
      pickerInput(
         NS(id, "year"),
         "Year",
         choices = all_year_options,
         selected = all_year_options,
         multiple = TRUE,
         options = pickerOptions(
           actionsBox = TRUE,
           selectedTextFormat = "count > 3"
         )
      ),
      pickerInput(
         NS(id, "condition"),
         "Condition",
         choices = all_condition_options,
         selected = all_condition_options,
         multiple = TRUE,
         options = pickerOptions(
           actionsBox = TRUE,
           selectedTextFormat = "count > 3"
         )
      ),
      pickerInput(
         NS(id, "background"),
         "Ancestor Strain",
         choices = all_background_options,
         selected = all_background_options,
         multiple = TRUE,
         options = pickerOptions(
           actionsBox = TRUE,
           selectedTextFormat = "count > 3"
         )
      ),
      pickerInput(
         NS(id, "sample"),
         "Sample",
         choices = all_sample_options,
         selected = all_sample_options,
         multiple = TRUE,
         options = pickerOptions(
           actionsBox = TRUE,
           selectedTextFormat = "count > 3"
         )
      ),
      ns = NS(id)
    ),
    # conditionalPanel(
    #   condition = "input.fhccSepView || output.selectedCumulView",
    #   ),
  )
}

region2gene_name <- function(gene_region, gene_info) {
  gene_info[gene_info$REGION == gene_region, "GENE"][1]
}

upload_vcf_data <- function(vcf_file, instructor, year, mut_backend, gene_info, session) {
      data <- read.csv(vcf_file$datapath)
      required_columns <- c(
        "CHROM", "POS", "REF", "ALT", "ANNOTATION", "REGION", "PROTEIN", "background", "condition", "sample"
      )
      if (!is.null(file) && all(required_columns %in% colnames(data))) {
        # Add instructor and year columns
        data$instructor <- rep(instructor, nrow(data))
        data$year <- rep(year, nrow(data))
        data$ID <- rep("NA", nrow(data))
        data$QUAL <- rep("NA", nrow(data))
        data$FILTER <- rep("NA", nrow(data))
        data$INFO <- rep("NA", nrow(data))
        data$seq_file <- rep("NA", nrow(data))

        if (!('GENE' %in% data)) {
          data$GENE <- sapply(data$REGION, function(region_name) {region2gene_name(region_name, gene_info)})
        }
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
        removeModal()
        sendSweetAlert(
          session = session,
          title = "Success",
          text = "VCF data uploaded successfully",
          type = "success"
        )
        return(df)
      } else {
        removeModal()
        sendSweetAlert(
          session = session,
          title = "Data Upload Error",
          text = "Some required columns are missing in your uploaded VCF file.",
          type = "error"
        )
        return(mut_backend)
      }
    }

selection_panel_server <- function(id, filtered_data, mutation_data, mut_backend, gene_info) {
  moduleServer(id, function(input, output, session) {
    mutation_data(mut_backend)
    observeEvent(input$uploadData, {
      showModal(modalDialog(
        title = "Upload VCF data",
        footer = modalButton("Cancel"),
        fileInput(NS(id, "datafile"), "Optional: Upload additional CSV File", accept = ".csv"),
        textInput(NS(id, "inputted_instructor"), "Who is your Instructor"),
      textInput(NS(id, "inputted_year"), "What is the current year"),
      actionButton(NS(id, "submit_teach_year"), "Submit Teacher and Year")
      ))
    })
    # Initialize a reactive variable for the dataframe
    # Function to read and append the uploaded data to the cumulative dataframe
    observeEvent(input$submit_teach_year, {
      mutation_data(upload_vcf_data(input$datafile, input$inputted_instructor, input$inputted_year, mutation_data(), gene_info, session))
      print(mutation_data()$instructor %>% unique)
    })

    # Display settings
    observe({
      if (!is.null(input$View)) {
        if (input$View == classroom_dropdown_text) {
          shinyjs::disable("fhccSepView")
          shinyjs::enable("classView")
        } else if (input$View == fhcc_sep_yevo_module_dropdown_text) {
          shinyjs::enable("fhccSepView")
          shinyjs::disable("classView")
        }
      } else {
        shinyjs::disable("classView")
        shinyjs::disable("fhccSepView")
      }
    })

    # Handling behaviors for button selections
    observe({
      options <- as.character(mutation_data() %>%
                                pull(background))
      # if (input$condition != "All Selected") {
      # options <- as.character(mutation_data() %>%
      #                           filter(condition == input$condition) %>%
      #                           pull(background))
      # }
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

    # observe({
    #   updateSelectInput(session, "instructor", choices = c(
    #     "All Selected",
    #     unique(mutation_data()$instructor)
    #   ))
    # })

    # observe({
    #   # Make sure the data exists before we do anything
    #   req(mutation_data())
    #
    #   # Figure out which years to show
    #   if (input$instructor == "All Selected") {
    #     years <- sort(unique(as.character(mutation_data()$year)))
    #   } else {
    #     years <- sort(unique(as.character(
    #       mutation_data()$year[
    #         mutation_data()$instructor == input$instructor
    #       ]
    #     )))
    #   }
    #
    # observe({
    #   updateSelectInput(session, "instructor", choices = c(
    #     "All Selected",
    #     unique(mutation_data()$instructor)
    #   ))
    # })

    # observe({
    #   updateSelectInput(session, "sample", choices = c(
    #     "All Selected",
    #     as.character(mutation_data() %>%
    #                    filter(instructor == input$instructor) %>%
    #                    filter(year == input$year) %>%
    #                    pull(sample))
    #   ))
    # })


    # download button functionality
    # output$downloadBtn <- downloadHandler(
    #   filename = function() {
    #     # Set the filename for the downloaded file
    #     if (is.null(input$View) ||
    #       (input$instructor == "All Selected" &&
    #         input$condition == "All Selected")) {
    #       "master_table.csv"
    #     } else {
    #       if (input$View == classroom_dropdown_text) {
    #         selected_instructor <- input$instructor
    #         selected_year <- input$year
    #         selected_sample <- input$sample
    #         if (selected_sample != "All Selected") {
    #           paste0(
    #             selected_instructor, "_", selected_year, "_",
    #             selected_sample, ".csv"
    #           )
    #         } else if (selected_year != "All Selected") {
    #           paste0(selected_instructor, "_", selected_year, ".csv")
    #         } else if (selected_instructor != "All Selected") {
    #           paste0(selected_instructor, ".csv")
    #         }
    #       } else {
    #         if (input$View == fhcc_sep_yevo_module_dropdown_text) {
    #           if (input$background != "All Selected") {
    #             paste0(input$condition, "_", input$background, ".csv")
    #           } else {
    #             paste0(input$condition, ".csv")
    #           }
    #         }
    #       }
    #     }
    #   },
    #   content = function(file) {
    #     # Write the data to a CSV file
    #     write.csv(filtered_data(), file)
    #   }
    # )

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
        value <- (input$View == classroom_dropdown_text)
      }
    })
    outputOptions(output, "selectedClassView", suspendWhenHidden = FALSE)

    output$selectedCumulView <- reactive({
      if (!is.null(input$View)) {
        value <- (input$View == fhcc_sep_yevo_module_dropdown_text)
      }
    })
    outputOptions(output, "selectedCumulView", suspendWhenHidden = FALSE)

    # filtering dataframe based on menu selection, most things from here on out
    # should be based on filtered_data()
    new_filtered_data <- reactive({
      data <- mutation_data()
      # filtering based on selections if NOT all selected
      # if (input$instructor != "All Selected") {
      #   data <- data %>% filter(instructor == input$instructor)
      # }
      # if (input$year != "All Selected") {
      #   data <- data %>% filter(year == input$year)
      # }
      # if (input$sample != "All Selected") {
      #   data <- data %>% filter(sample == input$sample)
      # }
      # if (input$condition != "All Selected") {
      #   data <- data %>% filter(condition == input$condition)
      # }
      # if (input$background != "All Selected") {
      #   data <- data %>% filter(background == input$background)
      # }
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