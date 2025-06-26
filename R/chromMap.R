chrom_map_ui <- function(id) {
  tabPanel(
    "Chromosome Map",
    plotlyOutput(NS(id, "chromPlot"), height = "600px"),
    verbatimTextOutput(NS(id, "info"))
  )
}

chrom_map_server <- function(
  id, final_gene, formatted_loading_message,
  chrom_info) {
  moduleServer(id, function(input, output, session) {
    # Displays Chromosome Map info; filtering by sample
    output$info <- renderText({
      return("Drag over variant tick mark to see details\n")
    })

    output$chromPlot <- renderPlotly({
      validate(
        need(final_gene()$START, formatted_loading_message)
      )

      # Preparing data
      final_gene_data <- final_gene()
      final_gene_data$CHROM <- factor(final_gene_data$CHROM,
                                      levels = levels(chrom_info$CHROM)
      )

      # Splitting data so that the colors don't get messed up with the single
      # counts
      final_gene_singletons <- final_gene_data %>% filter(Counts == 1)
      final_gene_multi_muts <- final_gene_data %>% filter(Counts >= 2)

      # Add numerical position for graphing chromosomes
      final_gene_singletons$chrom_num <- as.numeric(
        final_gene_singletons$CHROM
      ) - 1
      final_gene_multi_muts$chrom_num <- as.numeric(
        final_gene_multi_muts$CHROM
      ) - 1

      # Creating chromosome bars
      fig <- plot_ly(
        data = chrom_info,
        x = ~length,
        y = ~CHROM,
        type = "bar",
        orientation = "h",
        marker = list(color = "lightblue"),
        name = "Chromosomes",
        text = ~CHROM,
        textposition = "none",
        hoverinfo = "text",
        showlegend = FALSE
      )
      fig <- fig %>%
        layout(
          bargap = 0.3 # giving a little more space between chrom bars
        )


      # Shapes list for mutation markers
      shapes_list <- list()

      # Adding singleton mutation rectangles
      if (nrow(final_gene_singletons) > 0) {
        for (i in seq_len(nrow(final_gene_singletons))) {
          shapes_list[[length(shapes_list) + 1]] <- list(
            type = "rect",
            x0 = final_gene_singletons$START[i],
            x1 = final_gene_singletons$START[i] + 8000,
            y0 = final_gene_singletons$chrom_num[i] - 0.4,
            y1 = final_gene_singletons$chrom_num[i] + 0.4,
            fillcolor = "white",
            line = list(color = "black", width = 0.1),
            layer = "above"
          )
        }

        # Adding invisible points for singleton hover functionality
        fig <- fig %>%
          add_trace(
            data = final_gene_singletons,
            x = ~START + 4000,
            y = ~CHROM,
            type = "scatter",
            mode = "markers",
            marker = list(size = 10, opacity = 0, color = "white"),
            name = "Single Mutations",
            text = ~paste0(
              "Gene Name: ", GENE, "<br>Independent Mutations: ",
              Counts
            ),
            hoverinfo = "text"
          )
      }

      # Add multi-mutation rectangles
      if (nrow(final_gene_multi_muts) > 0) {
        # Create color scale
        max_count <- max(final_gene_multi_muts$Counts)
        # Create a vector to store colors for each multi mutation
        multi_mut_colors <- character(nrow(final_gene_multi_muts))

        for (i in seq_len(nrow(final_gene_multi_muts))) {
          # Calculate color based on count
          count_ratio <- (final_gene_multi_muts$Counts[i] - 2) / (max_count - 2)
          if (is.na(count_ratio) || count_ratio < 0) count_ratio <- 0

          # Linear interpolation between pink and dark red to match default
          # color interpolation on points
          r <- 255 - count_ratio * (255 - 139) # pink(255) to red4(139)
          g <- 192 - count_ratio * 192 # pink(192) to red4(0)
          b <- 203 - count_ratio * 203 # pink(203) to red4(0)

          hex_color <- rgb(r / 255, g / 255, b / 255)
          multi_mut_colors[i] <- hex_color

          shapes_list[[length(shapes_list) + 1]] <- list(
            type = "rect",
            x0 = final_gene_multi_muts$START[i],
            x1 = final_gene_multi_muts$START[i] + 8000,
            y0 = final_gene_multi_muts$chrom_num[i] - 0.4,
            y1 = final_gene_multi_muts$chrom_num[i] + 0.4,
            fillcolor = hex_color,
            line = list(color = "black", width = 0.1),
            layer = "above"
          )
        }
        final_gene_multi_muts$color <- multi_mut_colors

        # Adding invisible points for multi mutation hover functionality
        fig <- fig %>%
          add_trace(
            data = final_gene_multi_muts,
            x = ~START + 4000,
            y = ~CHROM,
            type = "scatter",
            mode = "markers",
            marker = list(
              color = ~Counts, # assuming you want to map 'Counts' to color
              colorscale = list(
                list(0, "#FFC0CB"),
                list(1, "#8B0000")
              ),
              size = 10,
              opacity = 0,
              colorbar = list(title = "Mutation Count")
            ),
            name = "Multiple Mutations",
            text = ~paste0(
              "Gene Name: ", GENE, "<br>Independent Mutations: ",
              Counts
            ),
            hoverinfo = "text"
          )
      }


      # Addding shapes to graph
      fig <- fig %>% layout(
        title = list(
          text = "Location of mutations along chromosomes",
          x = 0.5,
          xanchor = "center",
          y = 0.95,
          yanchor = "top",
          font = list(size = 20)
        ),
        xaxis = list(
          title = "Position along chromosome",
          showgrid = FALSE
        ),
        yaxis = list(
          title = "Chromosome",
          showgrid = FALSE,
          autorange = "reversed"
        ),
        shapes = shapes_list,
        plot_bgcolor = "rgba(0,0,0,0)",
        paper_bgcolor = "rgba(0,0,0,0)"
      )

      fig
    })
  })
}