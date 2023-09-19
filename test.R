
library(shiny)

ui <- fluidPage(
  titlePanel("Upload and Replace DataFrame"),
  sidebarLayout(
    sidebarPanel(
      fileInput("datafile", "Choose a CSV file"),
      actionButton("replaceBtn", "Replace DataFrame")
    ),
    mainPanel(
      tableOutput("data_table")
    )
  )
)
server <- function(input, output, session) {
  # Initialize a reactive variable for the dataframe
  uploaded_data <- reactiveVal(NULL)
  
  # Function to read and store the uploaded data as a dataframe
  observeEvent(input$datafile, {
    file <- input$datafile
    if (!is.null(file)) {
      df <- read.csv(file$datapath, sep = ",")
      uploaded_data(df)
    }
  })
  
  # Render the dataframe in the tableOutput
  output$data_table <- renderTable({
    uploaded_data()
  })
  
  # Replace the dataframe with the uploaded data
  observeEvent(input$replaceBtn, {
    if (!is.null(uploaded_data())) {
      # Perform any necessary data processing or validation here
      # For example, you can check if the uploaded data meets certain criteria.
      # If it does, replace the existing dataframe with the new one.
      
      # Replace the dataframe
      # You can also perform any additional data validation or cleaning here.
      uploaded_data(NULL)  # Clear the existing dataframe
    }
  })
}
shinyApp(ui, server)
