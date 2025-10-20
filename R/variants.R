variants_ui <- function(id) {
  tabPanel("Variant Pie Chart",
           plotlyOutput(NS(id, "varPieChart"),
                        height = "675px",
                        width = "100%"
           ), verbatimTextOutput(NS(id, "text")))
}

variants_server <- function(id, mutation_data, filtered_data, color_vector) {
  moduleServer(id, function(input, output, session) {
    output$varPieChart <- renderPlotly({

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
        ) %>%
        config(
          toImageButtonOptions = list(
            format = 'svg',
            filename = 'variant_pie_chart',
            height = 675,
            width = 800
          )
        )
    })
  })
}