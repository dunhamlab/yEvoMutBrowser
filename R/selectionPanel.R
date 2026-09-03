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
    div("", style = "height: 10px;"), # Create a 10px vertical space
    
    # create download button here
    downloadButton(NS(id, "downloadBtn"), "Download the following data"),
    div("", style = "height: 10px;"), # Create a 10px vertical space
    
    # A single combined set of filters -- no more "View By Class" vs "View By
    # Selection Condition" toggle. All five fields are always visible,
    # multi-select, and cross-filter each other symmetrically (see
    # selection_panel_server). "All Selected" is placeholder text shown
    # whenever a field is empty, rather than a real selectable choice.
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
    selectizeInput(NS(id, "condition"), "Condition", choices = "",
                   multiple = TRUE,
                   options = list(placeholder = "All Selected")),
    selectizeInput(NS(id, "background"), "Ancestor Strain", choices = "",
                   multiple = TRUE,
                   options = list(placeholder = "All Selected")),
  )
}

region2gene_name <- function(gene_region, gene_info) {
  gene_info[gene_info$REGION == gene_region, "GENE"][1]
}

# Treats a multi-select input as "unfiltered" if nothing is picked yet.
# ("All Selected" is placeholder text rather than a real choice, so "nothing
# picked" IS the all-selected state.)
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
    
    # ---- Symmetric (order-independent) cross-filtering --------------------
    # All five fields (instructor, year, sample, condition, background) are
    # cross-filtering multi-selects. Each field's choices come from
    # mutation_data() filtered by whatever is currently selected in the
    # OTHER four fields (isolated from its own value so it doesn't
    # re-trigger itself), unioned with its own current selection so an
    # existing pick is never silently cleared -- filtering only narrows
    # what NEW options are offered, never what's already chosen. This lets
    # the user start from any field -- sample, condition, year, whatever --
    # and have the rest narrow to match.
    #
    # Built as one factory instead of five near-identical blocks, since the
    # input id for each field matches its column name in mutation_data().
    filter_fields <- c("instructor", "year", "sample", "condition", "background")
    
    # Cache the last choices pushed to each widget, so we skip redundant
    # updateSelectizeInput calls (this + isolate() above is what avoids the
    # widgets flashing).
    last_choices <- reactiveValues()
    for (f in filter_fields) last_choices[[f]] <- character(0)
    
    lapply(filter_fields, function(field) {
      other_fields <- setdiff(filter_fields, field)
      observe({
        req(mutation_data())
        data <- mutation_data()
        for (other in other_fields) {
          sel <- input[[other]]
          if (!is_all_selected(sel)) {
            data <- data[data[[other]] %in% sel, , drop = FALSE]
          }
        }
        valid <- sort(unique(as.character(data[[field]])))
        current <- isolate(input[[field]])
        choices <- union(valid, current)
        
        if (!identical(last_choices[[field]], choices)) {
          updateSelectizeInput(session, field, choices = choices, selected = current)
          last_choices[[field]] <- choices
        }
      })
    })
    
    # download button functionality
    output$downloadBtn <- downloadHandler(
      filename = function() {
        # Build a compact label per field: "All" if unfiltered, otherwise
        # the picked values joined with "-". Fields left unfiltered are
        # dropped from the filename entirely to keep it readable.
        fmt <- function(x) {
          if (is_all_selected(x)) NULL else paste(x, collapse = "-")
        }
        parts <- Filter(Negate(is.null), list(
          fmt(input$instructor),
          fmt(input$year),
          fmt(input$sample),
          fmt(input$condition),
          fmt(input$background)
        ))
        if (length(parts) == 0) {
          "master_table.csv"
        } else {
          paste0(paste(unlist(parts), collapse = "_"), ".csv")
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
    
    # filtering dataframe based on menu selection, most things from here on out
    # should be based on filtered_data()
    new_filtered_data <- reactive({
      data <- mutation_data()
      # filtering based on selections if NOT all selected
      # (all five fields can each hold multiple picked values; an empty
      # selection means "All Selected" / no filter)
      if (!is_all_selected(input$instructor)) {
        data <- data %>% filter(instructor %in% input$instructor)
      }
      if (!is_all_selected(input$year)) {
        data <- data %>% filter(year %in% input$year)
      }
      if (!is_all_selected(input$sample)) {
        data <- data %>% filter(sample %in% input$sample)
      }
      if (!is_all_selected(input$condition)) {
        data <- data %>% filter(condition %in% input$condition)
      }
      if (!is_all_selected(input$background)) {
        data <- data %>% filter(background %in% input$background)
      }
      data
    })
    
    # NOTE: with the View toggle removed, there's no longer a notion of an
    # "incomplete" filter state -- every field defaults to "All Selected"
    # and any combination of picks is valid, so this is always TRUE now. If
    # whatever consumes form_complete downstream (e.g. gating mutation
    # selection) actually needs a stricter condition -- like "at least one
    # filter must be chosen" -- let me know and I'll tighten this back up.
    form_complete <- reactive({
      TRUE
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