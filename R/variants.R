
variants_ui <- function(id) {
  tabPanel(
    "Variant Pie Chart",
    plotlyOutput(
      NS(id, "varPieChart"),
      height = "675px",
      width = "100%"
    ),
    verbatimTextOutput(NS(id, "text"))
  )
}


variants_server <- function(id, mutation_data, filtered_data, color_vector) {
  moduleServer(id, function(input, output, session) {
    
    output$varPieChart <- renderPlotly({
      
      # Prepare filtered data
      pie_data <- filtered_data() %>%
        count(ANNOTATION, name = "count") %>%
        mutate(percent = count / sum(count) * 100) %>%
        arrange(ANNOTATION)
      
      # Use the shared annotation-to-color mapping.
      # Existing annotations retain their assigned colors.
      pie_colors <- unname(get_annotation_colors(as.character(pie_data$ANNOTATION)))
      
      # Fallback for annotations that do not yet have a color.
      missing_annotations <- unique(
        as.character(pie_data$ANNOTATION)[is.na(pie_colors)]
      )
      
      if (length(missing_annotations) > 0) {
        warning(
          "No color assigned for annotations: ",
          paste(missing_annotations, collapse = ", "),
          ". Add them to ANNOTATION_COLORS in the configuration."
        )
      }
      
      # Create pie chart
      plot_ly(
        data = pie_data,
        labels = ~ANNOTATION,
        values = ~percent,
        type = "pie",
        text = ~paste(ANNOTATION, ": ", round(percent, 2), "%"),
        hoverinfo = "text",
        textinfo = "text",
        marker = list(colors = pie_colors)
      ) %>%
        layout(
          title = list(
            text = "Percentage of Variants by Type",
            x = 1,
            y = 0.95,
            xanchor = "right",
            yanchor = "top"
          ),
          showlegend = TRUE,
          margin = list(l = 75, r = 200, b = 150, t = 50),
          plot_bgcolor = "rgba(0,0,0,0)",
          paper_bgcolor = "rgba(0,0,0,0)",
          xaxis = list(
            showgrid = FALSE,
            zeroline = FALSE,
            showticklabels = FALSE
          ),
          yaxis = list(
            showgrid = FALSE,
            zeroline = FALSE,
            showticklabels = FALSE
          )
        ) %>%
        config(
          toImageButtonOptions = list(
            format = "svg",
            filename = "variant_pie_chart",
            height = 675,
            width = 800
          )
        )
    })
  })
}
