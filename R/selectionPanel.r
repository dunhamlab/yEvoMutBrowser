selectionPanel <- function(mut_backend) {
  sidebarPanel(
    width = 3,
    fileInput("datafile", "Optional: Upload additional CSV File",
      accept = ".csv"
    ),
    conditionalPanel(
      # Asks for instructor/year so we can modify user uploaded file so we can
      # combine it with our current data
      condition = "output.filesUploaded",
      textInput("inputted_instructor", "Who is your Instructor"),
      textInput("inputted_year", "What is the current year"),
      actionButton("submit_teach_year", "Submit Teacher and Year")
    ),
    radioButtons("View", "Select an option:",
      choices = c("View By Class", "View By Selection Condition"),
      selected = character(0)
    ),
    div("", style = "height: 10px;"), # Create a 10px vertical space

    # create download button here
    downloadButton("downloadBtn", "Download the following data"),
    div("", style = "height: 10px;"), # Create a 10px vertical space

    # Only shows on condition observeEvent
    conditionalPanel(
      # links condition to button via button key
      condition = "input.uploadData || input.classView ||
      output.selectedClassView",
      selectInput("instructor", "Instructor", choices = c("")),
      selectInput("year", "Year", choices = c("")),
      # PLEASE NOTE: "sample" is called "Sample Name" within ui to make it
      # easier to understand, but the official name in the df is sample
      selectInput("sample", "Sample Name", choices = c("")),
    ),
    conditionalPanel(
      condition = "input.cumulView || output.selectedCumulView",
      selectInput("condition", "Condition", choices = c(
        "All Selected",
        mut_backend %>% count(condition) %>% pull(condition)
      )),
      selectInput("background", "Ancestor Strain", choices = c("")),
    ),
  )
}