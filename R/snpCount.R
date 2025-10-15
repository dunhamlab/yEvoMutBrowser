snp_count_ui <- function(id) {
  tabPanel("SNP Counts", plotlyOutput(NS(id, "snpCountPlot")))
}

snp_count_server <- function(id, filtered_data, color_vector) {
  moduleServer(id, function(input, output, session) {
    output$snpCountPlot <- renderPlotly({
      # Categorizing data into appropriate categories for plotting
      categorized_data <- filtered_data() %>%
        mutate(transition = paste0(REF, " to ", ALT)) %>%
        mutate(transition = if_else(nchar(transition) > 6,
                                    "Indel", transition
        )) %>%
        mutate(mutation_type = case_when(
          transition %in% c(
            "A to G", "G to A", "C to T",
            "T to C"
          ) ~ "Transition",
          transition %in% c(
            "A to T", "T to A", "C to G", "G to C", "A to C",
            "C to A", "T to G", "G to T"
          ) ~ "Transversion",
          TRUE ~ "Indel"
        )) %>%
        mutate(combined_group = case_when(
          transition %in% c("A to G", "T to C") ~ "A to G",
          transition %in% c("C to T", "G to A") ~ "C to T",
          transition %in% c("A to T", "T to A") ~ "A to T",
          transition %in% c("C to G", "G to C") ~ "C to G",
          transition %in% c("A to C", "T to G") ~ "A to C",
          transition %in% c("G to T", "C to A") ~ "C to A",
          TRUE ~ "Indel"
        ))

      # Summarize the data by the new group variable
      summarized_data <- categorized_data %>%
        group_by(combined_group, mutation_type) %>%
        summarise(count = n(), .groups = "drop")

      color_map <- set_names(color_vector, c( "Transition", "Transversion", "Indel"))

      desired_order <- c(
        "Indel", "A to G", "C to T", "A to T", "C to G",
        "A to C", "C to A"
      )

      summarized_data$combined_group <- factor(summarized_data$combined_group,
                                               levels = desired_order
      )

      # Plotting the data with coloring by categories
      p <- ggplot(summarized_data, aes(
        x = (combined_group), y = count,
        fill = mutation_type
      )) +
        aes(text = paste(
          combined_group, "\nCount:", count, "\nMutation Type:",
          mutation_type
        )) +
        geom_bar(position = "dodge", stat = "identity") +
        theme_bw() +
        scale_fill_manual(values = color_map) +
        labs(fill = "Mutation Type") +
        theme(
          legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          # Adjust angle for better visibility
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
        ) +
        ggtitle("Single Nucleotide Changes") +
        xlab("Mutation Type and SNP Call")

      p <- ggplotly(p, tooltip = "text") %>%
        config(
          toImageButtonOptions = list(
            format = 'svg',
            filename = 'snp_count_plot',
            height = 500,
            width = 800
          )
        )
    })
  })
}