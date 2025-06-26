variants_ui <- function(id) {
  tabPanel("Variant Pie Chart",
           plotlyOutput(NS(id, "varPieChart"),
                        height = "675px",
                        width = "100%"
           ), verbatimTextOutput(NS(id, "text")))
}

variants_server <- function(id, mutation_data, filtered_data) {
  moduleServer(id, function(input, output, session) {
    output$varPieChart <- renderPlotly({
      color_vector <- c(
        "#9edae5", "#17becf", "#dbdb8d", "#bcbd22", "#c7c7c7",
        "#e377c2", "#7f7f7f", "#f7b6d2", "#c49c94", "#8c564b",
        "#c5b0d5", "#9467bd", "#ff9896", "#d62728", "#98df8a",
        "#2ca02c", "#ffbb78", "#ff7f0e", "#aec7e8", "#1f77b4"
      )

      # Get all unique annotations in a fixed order
      all_unique_anno <- mutation_data() %>%
        distinct(ANNOTATION) %>%
        arrange(ANNOTATION) %>%
        pull(ANNOTATION)

      # Create named color map for Plotly marker
      color_map <- setNames(
        color_vector[seq_along(all_unique_anno)],
        all_unique_anno
      )

      # Prepare filtered data
      pie_data <- filtered_data() %>%
        count(ANNOTATION, name = "count") %>%
        mutate(percent = count / sum(count) * 100) %>%
        arrange(ANNOTATION)

      # Get color for each annotation in filtered data
      pie_colors <- unname(color_map[pie_data$ANNOTATION])

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
            text = "Percentage of Variants by Type", x = 1,
            y = 0.95, xanchor = "right", yanchor = "top"
          ),
          showlegend = TRUE,
          margin = list(l = 75, r = 200, b = 150, t = 50),
          plot_bgcolor = "rgba(0,0,0,0)",
          paper_bgcolor = "rgba(0,0,0,0)",
          xaxis = list(
            showgrid = FALSE, zeroline = FALSE,
            showticklabels = FALSE
          ),
          yaxis = list(
            showgrid = FALSE, zeroline = FALSE,
            showticklabels = FALSE
          )
        )
    })
  })
}