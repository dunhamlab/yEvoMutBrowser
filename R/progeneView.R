gene_pro_view_ui <- function(id) {
  ns <- NS(id)
  prefix <- ns("")

  tabPanel(
    "Gene View2", div("", style = "height: 10px;"),
    plotlyOutput(NS(id, "geneViewPlot"), width = "600px"), verbatimTextOutput(NS(id, "gene")),
    selectInput(NS(id, "geneSelectDropDown"), "Gene", choices = NULL),

    tags$div(
      id = "molstar-parent",
      style = "position: relative; width: 600px; height: 600px;",
      tags$canvas(
        id    = "molstar-canvas",
        style = "position: absolute; top: 0; left: 0; width: 100%; height: 100%;"
      ),
    )
    ,

    tags$script(HTML(sprintf("
      window.MY_MODULE_NS = '%s';
    ", prefix))),

    # ⬇️  place the script LAST so Shiny is ready
    tags$script(src = "static/molstar-custom.js"),
    verbatimTextOutput(NS(id, "resiinfo")),
    # actionButton(NS(id, 'mutations'), 'Mutations'),
    # plotlyOutput(NS(id, "mutplot"), width = "600px", height = "150px"),
    
    tags$div(class="domain_div",
      tags$div(class="mut_button_div",
    actionButton(NS(id, 'mutations'), 'Mutations'),
      ),

      tags$div(class="mut_plot_div",
    plotlyOutput(NS(id, "mutplot"), height = "140px"),
      ),
  
    ),

    
    tags$div(class="domain_div",
      tags$div(class="domain_button_div",
        actionButton(NS(id, 'domain'), 'Pfam Domains'),
      ),

      tags$div(class="domainplot_div",
        plotlyOutput(NS(id, "domainplot"), height = "65px"),
      ),
    
    ),
    
    includeCSS("R/www/styling.css"),

    tags$div(class = "container",
      # Inner div with class "header"
      tags$div(class = "protein-cont",
        tags$div(class = "my-verbatim-text",  # <-- your custom class
          textOutput(NS(id, "content"))
        )
      ),
      # Content section
      tags$div(class = "content_div",
        tags$div(class = 'url_div',
          tags$p("Protein: ", uiOutput(NS(id, "url"), inline = TRUE)),
        ),      
        tags$div(class = "description_div",
          tags$p("Description: ", textOutput(NS(id, "description"), inline = TRUE))
        ),

        tags$div(class = "pathway_div",
          tags$p("Pathways: ", uiOutput(NS(id, "pathway"), inline = TRUE))
        ),
        )
),

    # uiOutput(NS(id, "url"))
   )
}

gene_pro_view_server <- function(id, total_spaces, filtered_data, genes_info, link, gene_info_link_function, color_vector) {
  moduleServer(id, function(input, output, session) {


    mut_selected <- reactiveVal(FALSE)
    dom_selected <- reactiveVal(FALSE)
    #gene view draopdown menu
    observe({
      choices <- filtered_data() %>% pull(GENE)

      # Filter choices to include only those present in genes_info
      choices <- choices[choices %in% genes_info$GENE]

      choices <- sort(choices)

      updateSelectizeInput(session, "geneSelectDropDown",
                           choices = as.character(choices),
                           server = TRUE,
                           options = list(maxOptions = length(choices))
      )
    })

    gene <- reactive(input$geneSelectDropDown)


    cur_gene <- reactive({
      # merge & filter gene info and mutation info and filter for selected gene
      fd <- filtered_data() %>%
        merge(genes_info, by = intersect(names(filtered_data()), names(genes_info))) %>%
        arrange(GENE)

      filter(fd, GENE == gene())
    })

    # cg <- cur_gene()

    output$content <- renderText("Protein Info Card")

    output$description <- renderText(paste0(first(cur_gene()$DESCRIPTION)))

    all_annotations <- c(
      "missense", "nonsense", "5'-upstream",
      "indel-frameshift", "indel-inframe", "synonymous")

    letter_aa <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H",
      "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
    )

    three_letter_aa <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln",
                         "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe",
                         "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")

    names(three_letter_aa) <- letter_aa

    pfam <- read.csv(ORGANISM_PFAM_DOMAIN_INFO_PATH)
    pathway_info <- read.csv(ORGANISM_PATHWAY_INFO_PATH)

    output$pathway <- renderUI({
      validate(
        need(gene(), "LOADING"),
        need(filtered_data(), "Loading data...")
      )

      sgd <- first(cur_gene()$SGDID)
      row <- pathway_info[pathway_info$SGDID == sgd, ]
      BASE_URL <- "https://pathway.yeastgenome.org/YEAST/new-image?type=PATHWAY&object="

      if (nrow(row) == 0) {
                           tags$p("N/A")} else {

        links <- lapply(seq_len(nrow(row)), function(i) {
          url <- paste0(BASE_URL, row$PathwayID[i])

          tags$div(
            tags$a(href = url, target = "_blank", row$Pathway[i])
          )
        })

    # return a container with all links
    do.call(tags$div, links)
      }


  #     ifelse(nrow(row) == 0,
  #     paste0("N/Aa"),
  #     paste0(row$Pathway)
  #  )

    })

    genedatatable <- function(cur_gene) {
      pattern <- "(?<=\\d)([A-Za-z]|\\*|indel)$|([A-Za-z]|\\*)$"

      count_proteins <- cur_gene %>%
        mutate(indel = nchar(ALT) - nchar(REF)) %>% # Calculate indel difference
        group_by(POS, PROTEIN) %>%
        summarize(
          GENE = first(GENE),
          PROTEIN = first(PROTEIN),
          ANNOTATION = first(ANNOTATION),
          COUNTS = n(),
          Letter1 = substr(PROTEIN, 1, 1), # Extract the first character
          # Extract the numbers
          Numbers = as.numeric(str_extract(PROTEIN, "[0-9]+")),
          Letter2 = str_extract(PROTEIN, pattern),
          indel = first(indel)
        ) %>%
        ungroup()

      count_proteins_same <- count_proteins %>%
        group_by(Numbers) %>%
        summarize(
          GENE = first(GENE),
          PROTEIN = list(PROTEIN),
          ANNOTATION = first(ANNOTATION),
          Counts_diff_mutation = list(COUNTS),
          Counts_tot = sum(COUNTS),
          indel = first(indel)
        ) %>%
        ungroup()


      count_proteins_same <- count_proteins_same %>%
        mutate(
          combined = map2_chr(
            PROTEIN, Counts_diff_mutation,
            function(prot, counts) {
              prot_list <- str_split(prot, ", ") %>% unlist()
              counts_list <- str_split(counts, ", ") %>% unlist()

              combined_strings <- map2_chr(
                prot_list, counts_list,
                function(p, c) {
                  letter1 <- substr(p, 1, 1)
                  # numbers <- str_extract(p, "[0-9]+") %>% as.numeric()
                  letter2 <- str_extract(p, pattern)
                  paste("Count", three_letter_aa[letter1], "->", three_letter_aa[letter2], ":", c, "\n",
                    sep = " "
                  )
                }
              )
              paste(combined_strings, collapse = "")
            }
          )
        ) %>%
        mutate(
        PROTEIN = as.character(PROTEIN),
        # Extract the first character Amino Acid Wild Type
        AA_WT = substr(PROTEIN, 1, 1),
        AA_POS = if_else(ANNOTATION == "5'-upstream", -15,
            # Extract Amino Acid Position
            as.numeric(str_extract(PROTEIN, "[0-9]+"))
        ),
        # Amino Acid Mutation
        AA_M = substr(PROTEIN, nchar(PROTEIN), nchar(PROTEIN)),
        ANNOTATION = factor(ANNOTATION, levels = all_annotations)
        )


      count_proteins_same
    }


    # Learn about Gene button within gene viewer
    if(link != "NONE") {
    output$url <- renderUI({
      url <- a(gene(),
               href = gene_info_link_function(genes_info, input$geneSelectDropDown),
                target = "_blank"
      )
      url
    })
    }

    output$resiinfo <- renderText({
      req(input$resi_aa, input$resi_num)
      new <- genedatatable(cur_gene())
      residue_num <- input$resi_num
      if (residue_num %in% new$Numbers) {
        row <- new[new$Numbers == residue_num, ]

        ifelse(
          grepl("indel", row$ANNOTATION), # Check if ANNOTATION contains "indel"
          paste0(
              "Indel\n", abs(row$indel), " base ", ifelse(row$indel > 0,
                                                      "insertion", "deletion"
              ), "\nCount: ", row$Counts_diff_mutation,
              "\nPosition: ", row$AA_POS
          ), # Text for indel annotations
          ifelse(
            is.na(row$PROTEIN),
            paste0(
              row$ANNOTATION, "Count: ", row$Counts_diff_mutation,
              "Position: ", -abs(unique(row$POS) -
                                   unique(row$START))
            ),
            paste0(row$combined, "Position: ", row$AA_POS)
          )
        )
      }
      else {
      paste0("Position: ", input$resi_num, " Amino Acid: ", input$resi_aa)

      }
    })


    observeEvent(input$geneSelectDropDown, {
      gene_name <- input$geneSelectDropDown
      mut_selected(FALSE)
      dom_selected(FALSE)
      uniprotid <- genes_info$UniprotID[genes_info$GENE == gene_name]
      session$sendCustomMessage("initMolstar", uniprotid)
  })


    rect_data <- function(filtered_rows) {
      cg <- cur_gene()
      pd <- filtered_rows
      if (nrow(pd) == 0) {
        return(tibble(xmin=double(), xmax=double(),
                      ymin=double(), ymax=double(),
                      text=character(), id=character()))
      }
      rects <- tibble(
        xmin = pd$Start,
        xmax = pd$End,
        ymin = 0,
        ymax = 1,
        text = paste0(pd$Start, "-", pd$End, "\n", pd$Description)
      ) 
        
      rects$id <- paste0(first(cg$GENE), "_", seq_len(nrow(rects)), "_", rainbow(nrow(rects)))
      rects$fill_color <- rainbow(nrow(rects))
      rects

    }

    annotation_colors <- set_names(color_vector, all_annotations)

    new_theme_empty <- theme_bw()
    new_theme_empty$line <- element_blank()
    new_theme_empty$rect <- element_blank()
    new_theme_empty$strip.text <- element_blank()
    new_theme_empty$axis.text.y <- element_blank()
    new_theme_empty$axis.title <- element_blank()
    new_theme_empty$plot.margin <- structure(c(0, 0, -1, -1),
                                            unit = "lines",
                                            valid.unit = 3L,
                                            class = "unit")

    hover_text <- function(data, new_line) {
        ifelse(
          grepl("indel", data$ANNOTATION), # Check if ANNOTATION contains "indel"
          paste0(
            "Indel\n", abs(data$indel), " base ", ifelse(data$indel > 0,
                                                    "insertion", "deletion"
            ), "\nCount: ", data$Counts_diff_mutation,
            "Position: ", data$AA_POS
          ), # Text for indel annotations
          ifelse(
            is.na(data$PROTEIN),
            paste0(
              ANNOTATION, "\nCount: ", data$Counts_diff_mutation,
              "Position: ", -abs(unique(data$POS) -
                                     unique(data$START))
            ),
              paste0(data$combined, "Position: ", data$AA_POS)
            )
          )
      }

    
    output$geneViewPlot <- renderPlotly({
      ranges <- reactiveValues(x = NULL, y = NULL)

      # Showing please select gene message
      # Construct the string with spaces on each side of the loading message
      geneview_message <- "Please select a gene below"
      # Calculate the number of empty spaces needed on each side
      gene_message_length <- nchar(geneview_message)
      gene_spaces_on_each_side <- floor(
        (total_spaces - gene_message_length) / 2
      )
      select_gene_message <- sprintf(
        "%s%s%s%s",
        "\n\n\n\n\n\n\n\n",
        paste(rep(" ", gene_spaces_on_each_side), collapse = ""),
        geneview_message,
        paste(rep(" ", gene_spaces_on_each_side), collapse = "")
      )

      # Render gene view plot

      validate(
        need(gene(), select_gene_message),
        need(filtered_data(), "Loading data...")
      )

      cur_gene <- cur_gene()
      count_proteins_same <- genedatatable(cur_gene)
      # Pattern for extracting the second part of protein


      annotation_colors <- set_names(color_vector, all_annotations)

      # Generating ranges for the plot size
      xmax <- genes_info %>%
        filter(GENE == gene()) %>%
        pull(PROTEIN_LENGTH) %>%
        unique() %>%
        as.numeric()
      ranges$x <- c(-50, xmax + 50)

      ranges$y <- c(0, max(count_proteins_same$Counts_tot) +
                      ifelse(max(count_proteins_same$Counts_tot) < 5, 4, 1))

      p <- count_proteins_same %>%
        ggplot(
          aes(
            x = AA_POS,
            y = Counts_tot,
          )
        ) +
        # dummy geom_point to get legend showing all_annotations regardless of
        # what data is displayed
        geom_point(data = data.frame(
          PROTEIN = rep(NA_character_, length(all_annotations)),
          AA_WT = rep(NA_character_, length(all_annotations)),
          AA_POS = rep(NA_real_, length(all_annotations)),
          AA_M = rep(NA_character_, length(all_annotations)),
          ANNOTATION = factor(all_annotations, levels = all_annotations),
          Counts_diff_mutation = I(rep(
            list(numeric(0)),
            length(all_annotations)
          )),
          Counts_tot = rep(NA_integer_, length(all_annotations)),
          combined = rep(NA_character_, length(all_annotations))
        ), aes(y = 0, color = ANNOTATION, text = NULL), size = 2) +
        geom_hline(yintercept = 0, linetype = 2, alpha = .2, aes(text = NULL)) +
        geom_segment(aes(x = 0, xend = xmax, y = 0, yend = 0, text = NULL),
          size = 15, color = "cornflowerblue"
        ) +
        geom_segment(aes(
          x = AA_POS, xend = AA_POS, y = 0, yend = Counts_tot,
          text = NULL
        ), color = "gray") +
        geom_point(aes(
          x = AA_POS,
          y = Counts_tot,
          color = ANNOTATION,
          text = hover_text(count_proteins_same)
        ), size = 2) +
        ggtitle(as.character(gene())) +
        theme_classic() +
        theme(
          axis.text.x = element_text(),
          axis.title.y = element_text(),
          axis.text.y = element_text(),
          axis.ticks.y = element_line(),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(20, 20, 0, 20)
        ) +
        xlab("Amino acid position") +
        ylab("Mutation Count") +
        scale_y_continuous(breaks = function(x) {
          seq(floor(min(x)),
              ceiling(max(x)),
              by = 1
          )
        }) +
        scale_color_manual(values = annotation_colors) + # Manual color scale
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
        annotate("text",
                 x = 1, y = Inf,
                 label = "Drag over mutations to see more", hjust = 0, vjust = 2,
                 color = "black", size = 5
        ) +
        guides(color = guide_legend(title = "Annotation"))


      print(count_proteins_same)

      ggplotly(p, tooltip = "text")
    })

    output$mutplot <- renderPlotly({
      cur_gene <- cur_gene()
      count_proteins_same <- genedatatable(cur_gene)
      cg <- cur_gene()

            # Create a scatter plot with ggplot2
      gg <- ggplot(count_proteins_same, aes(x = count_proteins_same$Numbers, y = 0, key=count_proteins_same$ANNOTATION, color = count_proteins_same$ANNOTATION,
        text = hover_text(count_proteins_same)
      )) +
        scale_color_manual(values = annotation_colors) + 
        geom_point(size = 3.2, alpha = 0.3, shape=15 ) +
        # geom_point(color = "#5b91e2", size = 2.2) + # Add points with specified color and size
        labs(
             x = "X-axis Variable",             # Add x-axis label
             y = "Y-axis Variable") +  xlim(0, first(cur_gene()$PROTEIN_LENGTH)) +

        new_theme_empty + labs(x = "My X-axis Label") + guides(fill = FALSE) + theme(legend.position = "none")


      gg <- ggplotly(gg, tooltip = "text", source = "mutplot", dynamicTicks = TRUE) %>%
      layout(
              title = list(text = NULL),
      margin = list(l = 6, r = 6, t = 6, b = 20, pad = 0),

          xaxis = list(
            rangemode = "nonnegative",

            range = c(0, first(cg$PROTEIN_LENGTH)),   # <— force initial range [0, protein_length]
            autorange = FALSE         # <— prevent Plotly from auto-expanding it
          ) 
        ) %>%
      
        config(scrollZoom = TRUE)

      gg
    })


  observe({
    ed <- event_data("plotly_click", source="mutplot", priority = 'event')
    req(ed)
    print(ed)
    ann <- as.character(ed$key[[1]])
    hex <- annotation_colors[ann]
    print(hex)
    hex_val <- unname(hex)[1]
    positions <- c(ed$x)

    session$sendCustomMessage("highlightResidueWithSphere",
                              list(positions = positions,
                              colorHex=hex_val))
  })


  observeEvent(input$mutations, {
      mut_selected(!mut_selected())
      if (mut_selected()) {
        gene_info <- genedatatable(cur_gene())

        for (i in 1:nrow(gene_info)) {
          row <- gene_info[i, ]
          residue <- (row["Numbers"][[1]])
          mut_type <- row["ANNOTATION"][[1]]
          hex <- annotation_colors[mut_type]
          hex_val <- unname(hex)[1]
          session$sendCustomMessage("highlightResidueWithSphere",
                              list(positions = residue,
                              colorHex=hex_val
                              ))

        }
      }

      else {
        session$sendCustomMessage("clearSpheres", list())
      }
    })


    observeEvent(input$domain, {
      dom_selected(!dom_selected())
      if (dom_selected()) {
        # rects <- pfam_rectangles()
        cg <- cur_gene()
        pd <- pfam[pfam$UniparcID == first(cg$UniparcID), ]
        rects <- rect_data(pd)
        if (nrow(rects) != 0) {
        for (i in 1:nrow(rects)) {
          row <- rects[i, ]
          domain_start <- row['xmin'][[1]]
          domain_end <- row["xmax"][[1]]
          hex <- row["fill_color"][[1]]
          hex <- unname(hex)
          session$sendCustomMessage("highlightDomains",
         list(residueStart=domain_start, residueEnd=domain_end, colorHex=hex))

        }}
      }

    else {
      session$sendCustomMessage("clearPaint", list())
    }
    })


    plot_rect <- function(domain_db) {
      # rects <- pfam_rectangles()
      cg <- cur_gene()

      pd <- pfam[pfam$UniparcID == first(cg$UniparcID), ]
      rects <- rect_data(pd)

      if (nrow(rects) == 0) {
        return(plotly_empty(type="scatter", mode="markers") %>%
                 layout(title = "No Pfam domains"))
      }


      p <- ggplot(rects) +
        geom_rect(aes(xmin = xmin, xmax = xmax,
                    ymin = ymin, ymax = ymax,
                    text = text,
                    key = id,
                    fill = fill_color
                  ),

                  color = "black", 
                  alpha = 0.6) + scale_fill_identity() +
        # labs(title = paste("Pfam Domains for", gene()), x = "Residue", y = "") +
        theme_minimal() + scale_x_continuous(limits = c(0, first(cg$PROTEIN_LENGTH))) + new_theme_empty +
        guides(fill = FALSE)


      gg <- ggplotly(p, tooltip = "text", dynamicTicks = TRUE, source = "domainplot") %>%
        layout(
          title = list(text = NULL),
          margin = list(l = 6, r = 6, t = 6, b = 20, pad = 0),

          xaxis = list(
            rangemode = "nonnegative",

            range = c(0, first(cg$PROTEIN_LENGTH)),    # <— force initial range [0, protein_length]
            autorange = FALSE         # <— prevent Plotly from auto-expanding it
          ) 
        ) %>% config(scrollZoom = TRUE)

      gg

    }


    output$domainplot <- renderPlotly({
        plot_rect(pfam)      
    })

    observe({
      # priortiy event will be triggered at every event/click
      ed <- event_data("plotly_click", source = "domainplot", priority = 'event')
      req(ed$key)
      print(ed)
      cg <- cur_gene()
      print(ed$key)
      str <- ((ed$key)[[1]])
      split_vector <- unlist(strsplit(str, "_"))
      hex <- tail(split_vector, 1)

      pd <- pfam[pfam$UniparcID == first(cg$UniparcID), ]
      rects <- rect_data(pd)
      domain <- rects[rects$id == ed$key, ]
      domain_start <- domain$xmin
      domain_end <- domain$xmax

      
      session$sendCustomMessage("highlightDomains",
      list(residueStart=domain_start, residueEnd=domain_end, colorHex=hex))
    })

  })
}